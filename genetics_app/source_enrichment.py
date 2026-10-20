"""Prepare public rsID placements independently of any personal DNA file.

Only variants selected from gaps in the public SNPedia catalogue and public
rules are requested. Analysis never calls this API or sends marker lists.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor,as_completed
from datetime import datetime,timezone
import json
from pathlib import Path
import re
import sqlite3
import time
import urllib.request
import urllib.error
from .references import digest

API='https://grch37.rest.ensembl.org/variation/homo_sapiens'
SELECTION_VERSION=1


def validate_cache(state,source,criteria,chain_path,complete=True):
    expected=dict(schema=1,kind='Ensembl public placements',build='37',source=API,
        database_sha256=digest(source),criteria_sha256=digest(criteria),
        chain_sha256=digest(chain_path) if chain_path else '',selection_version=SELECTION_VERSION)
    if any(state.get(k)!=v for k,v in expected.items()):raise ValueError('Incompatible public placement cache identity.')
    requested=state.get('requested',[]);records=state.get('records',{});missing=state.get('missing',[])
    if not isinstance(records,dict) or len(requested)!=len(set(requested)) or len(missing)!=len(set(missing)) or len(requested)>25000:raise ValueError('Invalid public placement inventory.')
    if not all(re.fullmatch(r'rs\d+',k) for k in requested):raise ValueError('Invalid public variant identifier.')
    known=set(records)|set(missing)
    if set(records)&set(missing) or not known<=set(requested):raise ValueError('Inconsistent public placement inventory.')
    if complete and (state.get('complete') is not True or known!=set(requested)):raise ValueError('Incomplete public placement cache.')
    for key,r in records.items():
        if not isinstance(r,dict) or key not in {r.get('name'),*r.get('synonyms',[])}:raise ValueError('Public placement identity is not corroborated.')
    return state


def response_records(data,batch):
    if not isinstance(data,dict) or 'error' in data:raise ValueError('Unexpected public placement response.')
    wanted=set(batch);indexed={}
    for key,record in data.items():
        if not isinstance(record,dict):raise ValueError('Malformed public variant record.')
        identities={key,record.get('name')}|set(record.get('synonyms',[]))
        related=identities&wanted
        if not related:raise ValueError('Response includes an unrelated variant.')
        for requested in related:
            if requested in indexed and indexed[requested]!=record:raise ValueError('Conflicting canonical aliases in public response.')
            indexed[requested]=record
    return indexed


def candidates(source,criteria,chain_path):
    from .knowledge import parameters,Chain
    chain=Chain(chain_path) if chain_path else None
    db=sqlite3.connect(Path(source).resolve().as_uri()+'?mode=ro',uri=True)
    wanted=set();known=set()
    for key,text in db.execute('select rsid,content from snps'):
        known.add(key)
        if not re.fullmatch(r'rs\d+',key):continue
        f=parameters(text,{'rsnum','23andme snp'});assembly=f.get('assembly',f.get('genomebuild',''))
        build='38' if '38' in assembly else '37' if '37' in assembly else None
        pos=int(f.get('position','0')) if f.get('position','0').isdigit() else 0
        ch=f.get('chromosome','').removeprefix('chr')
        mapped=chain.map(ch,pos,reverse=True) if chain and build=='38' and ch and pos else None
        if not build or not pos or not ch or build=='37' or not f.get('stabilizedorientation',f.get('orientation')) or f.get('flip37','').lower()=='true' or f.get('flip38','').lower()=='true' or (build=='38' and (not mapped or mapped[1]=='-')):wanted.add(key)
    for key, in db.execute('select id from genotypes'):
        key=key.split('(')[0]
        if key not in known and re.fullmatch(r'rs\d+',key):wanted.add(key)
    db.close()
    rules=json.loads(Path(criteria).read_text())
    for page in rules['pages']:
        for rev in page.get('revisions',[]):wanted.update(re.findall(r'\brs\d+',rev['slots']['main']['content']))
    return sorted(wanted)


def download(source,criteria,chain_path,output,progress=print):
    source,output=Path(source),Path(output)
    wanted=candidates(source,criteria,chain_path)
    if not wanted or len(wanted)>25000:raise ValueError('Unexpected public reference inventory size.')
    identity=dict(database_sha256=digest(source),criteria_sha256=digest(criteria),chain_sha256=digest(chain_path) if chain_path else '',selection_version=SELECTION_VERSION)
    state=dict(schema=1,kind='Ensembl public placements',build='37',source=API,requested=wanted,records={},missing=[],complete=False,**identity)
    if output.exists():
        previous=json.loads(output.read_text())
        validate_cache(previous,source,criteria,chain_path,complete=False)
        if any(previous.get(k)!=v for k,v in identity.items()) or previous.get('requested')!=wanted:raise ValueError('Public placement checkpoint belongs to a different source selection.')
        state=previous
    pending=[k for k in wanted if k not in state['records'] and k not in state['missing']]
    def fetch(batch):
        request=urllib.request.Request(API,data=json.dumps({'ids':batch}).encode(),headers={'Content-Type':'application/json','Accept':'application/json','User-Agent':'DNAResearchReport/1.0 public-reference-preparation'})
        for attempt in range(3):
            try:
                with urllib.request.urlopen(request,timeout=45) as r:data=json.load(r)
                break
            except (urllib.error.HTTPError,urllib.error.URLError,TimeoutError) as e:
                if isinstance(e,urllib.error.HTTPError) and e.code not in {429,500,502,503,504}:raise
                if attempt==2:raise
                time.sleep(3*(attempt+1))
        return batch,response_records(data,batch)
    # Bounded parallel public requests; only this thread mutates the checkpoint.
    with ThreadPoolExecutor(max_workers=3) as pool:
        jobs=[pool.submit(fetch,pending[start:start+200]) for start in range(0,len(pending),200)]
        for job in as_completed(jobs):
            batch,data=job.result()
            for key in batch:
                if key in data:state['records'][key]=data[key]
                else:state['missing'].append(key)
            state['retrieved_at']=datetime.now(timezone.utc).isoformat();state['complete']=len(state['records'])+len(state['missing'])==len(wanted)
            output.parent.mkdir(parents=True,exist_ok=True);temp=output.with_suffix('.tmp');temp.write_text(json.dumps(state,sort_keys=True));temp.replace(output)
            progress(f"Prepared {len(state['records'])+len(state['missing'])}/{len(wanted)} public rsID records")
    return state


def unique_snv(record,build='37'):
    """Unique forward primary-chromosome SNV; alternate scaffolds are ignored."""
    if record.get('var_class')!='SNP':return None
    found=set()
    primary={str(i) for i in range(1,23)}|{'X','Y','MT'}
    for m in record.get('mappings',[]):
        if m.get('assembly_name')!='GRCh'+build or m.get('seq_region_name') not in primary:continue
        alleles=m.get('allele_string','').split('/')
        if m.get('strand')!=1 or not isinstance(m.get('start'),int) or m['start']<1 or m.get('start')!=m.get('end') or not all(len(a)==1 and a in 'ACGT' for a in alleles):return None
        found.add((m['seq_region_name'],m['start'],tuple(sorted(set(alleles)))))
    if len(found)!=1:return None
    ch,pos,alleles=next(iter(found))
    if not 2<=len(alleles)<=4:return None
    return dict(chrom=ch,pos=pos,alleles=list(alleles),build=build,canonical_id=record.get('name'),synonyms=[x for x in record.get('synonyms',[]) if re.fullmatch(r'rs\d+',x)],source=API)


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--database',type=Path,required=True);p.add_argument('--criteria',type=Path,required=True);p.add_argument('--chain',type=Path);p.add_argument('--out',type=Path,required=True);a=p.parse_args()
    if a.out.resolve().is_relative_to(Path(__file__).resolve().parents[1]):p.error('Keep reference caches outside the release folder.')
    download(a.database,a.criteria,a.chain,a.out)

def enrich_marker(marker,placement,genotype_alleles):
    """Attach a public placement and derive strand only from complete alleles.

    No observed genotype or saved report participates in this operation.
    Palindromic or incomplete allele sets cannot determine an orientation.
    """
    if not placement or (marker['chrom'] and marker['chrom']!=placement['chrom']):return marker
    if marker.get('indel') or any(len(a)!=1 or a not in 'ACGT' for a in marker.get('alleles',[])):return marker
    build=placement['build'];before=marker['positions'].get(build)
    marker['public_placement']=placement
    marker['chrom']=placement['chrom'];marker['positions'][build]=placement['pos']
    source=set(marker.get('alleles',[]))|set(genotype_alleles)
    marker['source_snv_alleles']=sorted(source)
    forward=set(placement['alleles']);comp={a.translate(str.maketrans('ACGT','TGCA')) for a in source}
    matches=[name for name,alleles in [('plus',source),('minus',comp)] if alleles<=forward]
    if len(source)==2 and source<=set('ACGT') and len(matches)==1:
        marker.setdefault('orientations_by_build',{})[build]=matches[0]
    marker['public_placement_audit']=dict(original_position=before,original_orientation=marker['orientation'],
        method='unique_Ensembl_primary_SNV_placement',orientation_method='non_palindromic_public_source_allele_equivalence' if len(matches)==1 and len(source)==2 else 'not_inferred',
        source_alleles=sorted(source),forward_alleles=placement['alleles'])
    return marker


if __name__=='__main__':main()
