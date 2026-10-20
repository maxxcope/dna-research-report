"""Versioned SNPedia knowledge import and coordinate-aware local interpretation."""
import bisect
from collections import Counter, defaultdict
from datetime import datetime, timezone
import gzip
import hashlib
import html
import json
from pathlib import Path
import re
import sqlite3
import tempfile
from .references import digest, connect_read
from .rules import evaluate_all, parse_rule, unresolved_markers, conditional_evidence
from .variants import public_indels, source_alleles, prepare_indel, resolve_indel, encoding_warning

LICENSE = 'CC BY-NC-SA 3.0 US'

def templates(text):
    """Yield outer templates; split parameters only outside nested templates/links."""
    starts=[]; depth=0; i=0
    while i<len(text)-1:
        pair=text[i:i+2]
        if pair=='{{':
            if depth==0:start=i+2
            depth+=1;i+=2;continue
        if pair=='}}' and depth:
            depth-=1
            if depth==0:
                body=text[start:i]; parts=[]; token=''; nested=links=0; j=0
                while j<len(body):
                    pair=body[j:j+2]
                    if pair in ('{{','[['):
                        if pair=='{{':nested+=1
                        else:links+=1
                        token+=pair;j+=2;continue
                    if pair in ('}}',']]'):
                        if pair=='}}':nested=max(0,nested-1)
                        else:links=max(0,links-1)
                        token+=pair;j+=2;continue
                    if body[j]=='|' and not nested and not links:parts.append(token.strip());token=''
                    else:token+=body[j]
                    j+=1
                parts.append(token.strip())
                yield parts[0].lower(),parts[1:]
            i+=2;continue
        i+=1

def parameters(text,name):
    for title,parts in templates(text):
        if title in name:
            return {key.strip().lower():value.strip() for p in parts if '=' in p for key,value in [p.split('=',1)]}
    return {}

def plain(text):
    # Metadata is rendered in labelled report fields rather than as wiki markup.
    spans=[];depth=0;start=0
    for match in re.finditer(r'\{\{|\}\}',text):
        if match.group()=='{{':
            if depth==0:start=match.start()
            depth+=1
        elif depth:
            depth-=1
            if depth==0:
                raw=text[start:match.end()];parts=next(templates(raw),('',[]));name,values=parts
                if name in {'rsnum','23andme snp','genotype','genoset','population diversity','clinvar','omim','on chip'}:replacement=''
                elif name.startswith('pmid'):
                    pmid=re.search(r'\d{4,9}',values[0] if values else '')
                    replacement=' [Publication PMID '+pmid.group()+'] ' if pmid else ''
                else:continue
                spans.append((start,match.end(),replacement))
    for start,end,replacement in reversed(spans):text=text[:start]+replacement+text[end:]
    text=re.sub(r'<!--.*?-->','',text,flags=re.S)
    text=re.sub(r'\[\[([^\]|]+)\|([^\]]+)\]\]',r'\2',text)
    text=re.sub(r'\[\[([^\]]+)\]\]',r'\1',text)
    text=re.sub(r'\[(https?://\S+)\s+([^\]]+)\]',r'\2 (\1)',text)
    # Preserve source prose and literal template text as escaped text, never HTML.
    text=re.sub(r'<[^>]*>','',text).replace("'''",'').replace("''",'')
    return html.unescape(text).strip()

def bibliography(text):
    return sorted(set(re.findall(r'(?:PMID\s*(?:Auto(?: GWAS)?)?\s*[|=]|PMID\s*[=:]?\s*)(\d{4,9})',text,re.I)))

def number(value,default=0):
    try:return float(value)
    except (ValueError,TypeError):return default

class Chain:
    """One-to-one base mapping through same-chromosome UCSC hg19 -> hg38 blocks."""
    def __init__(self,path):
        self.forward=defaultdict(list);self.reverse=defaultdict(list)
        with gzip.open(path,'rt') as stream:
            for line in stream:
                f=line.split()
                if not f:continue
                if f[0]=='chain':
                    _,score,tc,ts,tstrand,tstart,tend,qc,qs,qstrand,qstart,qend,cid=f
                    tp,qp,qs=int(tstart),int(qstart),int(qs)
                    use=tc==qc and tstrand=='+';ch=tc.removeprefix('chr');continue
                size=int(f[0])
                if use:
                    qlo=qp if qstrand=='+' else qs-qp-size
                    block=(tp,tp+size,qlo,qlo+size,qstrand)
                    self.forward[ch].append(block);self.reverse[ch].append(block)
                if len(f)==3:tp+=size+int(f[1]);qp+=size+int(f[2])
        for ch in self.forward:
            self.forward[ch].sort(key=lambda x:x[0]);self.reverse[ch].sort(key=lambda x:x[2])
        self.starts={False:{ch:[x[0] for x in a] for ch,a in self.forward.items()},True:{ch:[x[2] for x in a] for ch,a in self.reverse.items()}}
        self.prefix_ends={}
        for rev,blocks in [(False,self.forward),(True,self.reverse)]:
            self.prefix_ends[rev]={}
            for ch,arr in blocks.items():
                mx=0;ends=[]
                for b in arr:mx=max(mx,b[3] if rev else b[1]);ends.append(mx)
                self.prefix_ends[rev][ch]=ends
    def map(self,ch,pos,reverse=False):
        arr=(self.reverse if reverse else self.forward).get(ch,[]);p=pos-1;hits=set()
        i=bisect.bisect_right(self.starts[reverse].get(ch,[]),p)-1
        while i>=0 and self.prefix_ends[reverse][ch][i]>p:
            t0,t1,q0,q1,strand=arr[i];lo,hi=(q0,q1) if reverse else (t0,t1)
            if lo<=p<hi:
                off=p-lo
                out=(t0+off if strand=='+' else t1-1-off) if reverse else (q0+off if strand=='+' else q1-1-off)
                hits.add((out+1,strand))
            i-=1
        return next(iter(hits)) if len(hits)==1 else None

def population_rows(text):
    out=[]
    for name,parts in templates(text):
        if name!='population diversity':continue
        genos={};items=[]
        for part in parts:
            if '=' in part:
                k,v=part.split('=',1);genos[k.strip().lower()]=v.strip()
            else:items.append(part.strip())
        for i in range(0,len(items)-3,4):
            label,*values=items[i:i+4]
            try:values=[float(v) for v in values]
            except ValueError:continue
            if not label or any(v<0 or v>100 for v in values):continue
            out.append(dict(population=label,genotypes=[genos.get('geno'+str(j),'') for j in (1,2,3)],percentages=values,release=genos.get('hapmaprevision','unspecified')))
    return out

def public_snv_positions(path, wanted):
    """Unique rsID locations from a versioned ClinVar SNV catalogue, never a genome.

    Exclude IDs with any non-SNV or multiple primary placements. Do not apply
    VCF indel anchor positions to array probes.
    """
    db, meta = connect_read(path)
    if meta.get('kind') != 'ClinVar' or meta.get('build') not in {'37','38'}:
        db.close(); raise ValueError('Coordinate corroboration requires a built ClinVar catalogue.')
    loci=defaultdict(set); rejected=set()
    try:
        for ch,pos,ref,alt,rs in db.execute('select chrom,pos,ref,alt,rs from variants'):
            for token in re.split(r'[;,]',rs):
                key='rs'+token
                if key not in wanted:continue
                if len(ref)!=1 or ref not in 'ACGT' or any(len(a)!=1 or a not in 'ACGT.' for a in alt.split(',')):
                    rejected.add(key);continue
                loci[key].add((ch,pos))
    finally:db.close()
    return {k:next(iter(v)) for k,v in loci.items() if len(v)==1 and k not in rejected},meta['build']

def import_knowledge(source,criteria,chain_path,destination,coordinate_reference=None,placement_reference=None,progress=lambda _:None):
    source,criteria,destination=map(Path,(source,criteria,destination))
    if destination.resolve() in (source.resolve(),criteria.resolve()):raise ValueError('Reference destination must differ from inputs.')
    rules=json.loads(criteria.read_text())
    if not rules.get('complete'):raise ValueError('Complete public criteria inventory required.')
    if rules.get('database_sha256') and rules['database_sha256']!=digest(source):raise ValueError('Criteria inventory belongs to a different source snapshot.')
    rule_by_id={}
    for page in rules['pages']:
        if 'revisions' not in page:continue
        rev=page['revisions'][0];content=rev['slots']['main']['content']
        if hashlib.sha1(content.encode()).hexdigest()!=rev['sha1']:raise ValueError('Criteria checksum mismatch.')
        parse_rule(content)
        rule_by_id[page['title'].split('/')[0].lower()]=(content,rev['revid'],rev['timestamp'])
    chain=Chain(chain_path) if chain_path else None
    src=sqlite3.connect(source.resolve().as_uri()+'?mode=ro',uri=True)
    if src.execute('pragma quick_check').fetchone()[0]!='ok':raise ValueError('Reference integrity check failed.')
    progress('Corroborating public source coordinates and allele definitions')
    corroborated,reference_build=public_snv_positions(coordinate_reference,{r[0] for r in src.execute('select rsid from snps')}) if coordinate_reference else ({},None)
    indels,indel_build=public_indels(coordinate_reference,{r[0] for r in src.execute('select rsid from snps')})
    from .source_enrichment import unique_snv,enrich_marker,validate_cache
    from .curation import definitions,apply_definition
    curated,curation_hash=definitions();curated={r['id']:r for r in curated}
    placements={};genotype_alleles=defaultdict(set)
    for key, in src.execute('select id from genotypes'):
        match=re.fullmatch(r'((?:rs|i)\d+)\(([acgt]);([acgt])\)',key,re.I)
        if match:genotype_alleles[match[1].lower()].update([match[2].upper(),match[3].upper()])
    if placement_reference:
        extra=json.loads(Path(placement_reference).read_text())
        validate_cache(extra,source,criteria,chain_path)
        placements={k:p for k,r in extra['records'].items() if (p:=unique_snv(r))}
    destination.parent.mkdir(parents=True,exist_ok=True)
    with tempfile.TemporaryDirectory(dir=destination.parent) as tmp:
        db=sqlite3.connect(str(Path(tmp)/'knowledge.sqlite'))
        db.executescript('CREATE TABLE metadata(key TEXT PRIMARY KEY,value TEXT); CREATE TABLE markers(id TEXT PRIMARY KEY,data TEXT); CREATE TABLE genotypes(id TEXT PRIMARY KEY,data TEXT); CREATE TABLE genosets(id TEXT PRIMARY KEY,data TEXT);')
        counts=Counter()
        marker_total=src.execute('select count(*) from snps').fetchone()[0]
        for marker,text,retrieved in src.execute('select rsid,content,scraped_at from snps'):
            f=parameters(text,{'rsnum','23andme snp'});ch=f.get('chromosome','').removeprefix('chr');pos=int(f.get('position','0')) if f.get('position','').isdigit() else 0
            build='38' if '38' in f.get('assembly',f.get('genomebuild','')) else '37' if '37' in f.get('assembly',f.get('genomebuild','')) else None
            coords={};mapping_strand='+'
            if ch and pos and build:coords[build]=pos
            # Official GRCh37 and GRCh38 both use NC_012920.1 for MT.
            # Do not use UCSC hg19 chrM's different mitochondrial sequence.
            if ch=='MT' and pos and build:coords={'37':pos,'38':pos}
            if chain and ch and ch!='MT' and pos and build:
                mapped=chain.map(ch,pos,reverse=build=='38')
                if mapped:coords['37' if build=='38' else '38']=mapped[0];mapping_strand=mapped[1]
            coordinate_audit=None
            verified=corroborated.get(marker)
            if verified and (not ch or ch==verified[0]):
                verified_ch,verified_pos=verified
                if coords.get(reference_build)!=verified_pos:
                    coordinate_audit=dict(method='unique_ClinVar_SNV_rsID_placement',build=reference_build,original_positions=coords.copy(),original_chrom=ch)
                    ch=verified_ch;coords={reference_build:verified_pos};mapping_strand='+'
                    if ch=='MT':coords={'37':verified_pos,'38':verified_pos}
                    elif chain:
                        mapped=chain.map(ch,verified_pos,reverse=reference_build=='38')
                        if mapped:coords['37' if reference_build=='38' else '38']=mapped[0];mapping_strand=mapped[1]
            record=dict(id=marker,chrom=ch,positions=coords,mapping_strand=mapping_strand,orientation=f.get('stabilizedorientation',f.get('orientation','')).lower(),flip37=f.get('flip37','').lower()=='true',flip38=f.get('flip38','').lower()=='true',alias=('rs'+f['rsid']) if marker.startswith('i') and f.get('rsid','').isdigit() else None,gene=plain(f.get('gene',f.get('gene_s',''))),text=plain(text),pmids=bibliography(text),populations=population_rows(text),retrieved=retrieved,source='https://bots.snpedia.com/index.php/'+marker)
            record['coordinate_audit']=coordinate_audit
            record['alleles']=source_alleles(f)
            record['indel']=prepare_indel(record,indels.get(marker),indel_build) if indel_build else None
            record=enrich_marker(record,placements.get(marker),genotype_alleles.get(marker,set()))
            if marker in curated:record=apply_definition(record,text,genotype_alleles.get(marker,set()),curated[marker])
            db.execute('insert into markers values(?,?)',(marker,json.dumps(record)));counts['markers']+=1
            if counts['markers'] % 2500 == 0:progress(f"Source markers indexed: {counts['markers']:,}/{marker_total:,}")
        known={r[0] for r in src.execute('select rsid from snps')}
        for marker in sorted((set(genotype_alleles)&placements.keys())-known):
            record=dict(id=marker,chrom='',positions={},mapping_strand='+',orientation='',flip37=False,flip38=False,alias=None,gene='',text='',pmids=[],populations=[],retrieved=extra.get('retrieved_at',''),source='https://bots.snpedia.com/index.php/'+marker,alleles=[],indel=None,coordinate_audit=None)
            enrich_marker(record,placements[marker],genotype_alleles[marker])
            db.execute('insert into markers values(?,?)',(marker,json.dumps(record)));counts['markers']+=1
        for key,text,retrieved in src.execute('select id,content,scraped_at from genotypes'):
            f=parameters(text,{'genotype'})
            record=dict(id=key,summary=plain(f.get('summary','')),magnitude=number(f.get('magnitude')),repute=f.get('repute','Unknown'),text=plain(text),pmids=bibliography(text),retrieved=retrieved,source='https://bots.snpedia.com/index.php/'+key.split('(')[0]+'('+key.split('(',1)[1].upper() if '(' in key else 'https://bots.snpedia.com/index.php/'+key)
            db.execute('insert into genotypes values(?,?)',(key,json.dumps(record)));counts['genotypes']+=1
            if counts['genotypes'] % 10000 == 0:progress(f"Source genotype explanations indexed: {counts['genotypes']:,}")
        for key,text,retrieved in src.execute('select id,content,scraped_at from genosets'):
            f=parameters(text,{'genoset'});r=rule_by_id.get(key)
            record=dict(id=key,summary=plain(f.get('summary','')),magnitude=number(f.get('magnitude')),repute=f.get('repute','Unknown'),text=plain(text),pmids=bibliography(text),rule=r[0] if r else None,revision=r[1] if r else None,rule_date=r[2] if r else None,source='https://bots.snpedia.com/index.php/'+key,retrieved=retrieved)
            db.execute('insert into genosets values(?,?)',(key,json.dumps(record)));counts['genosets']+=1
        meta=dict(schema='1',kind='SNPedia',records=str(counts['genotypes']),release_date=rules.get('retrieved_at',''),source_sha256=digest(source),criteria_sha256=digest(criteria),chain_sha256=digest(chain_path) if chain_path else '',license=LICENSE,license_url='https://creativecommons.org/licenses/by-nc-sa/3.0/us/',source_url='https://bots.snpedia.com/',prepared_utc=datetime.now(timezone.utc).isoformat(),counts=json.dumps(counts))
        meta.update(importer_version='6',coordinate_reference_sha256=digest(coordinate_reference) if coordinate_reference else '',placement_reference_sha256=digest(placement_reference) if placement_reference else '',curated_metadata_sha256=curation_hash)
        db.executemany('insert into metadata values(?,?)',meta.items());db.commit();db.close();src.close();Path(tmp,'knowledge.sqlite').replace(destination)
    return meta

def interpret(genome,path,flagged=frozenset()):
    db,meta=connect_read(path)
    try:
        markers={r['id']:json.loads(r['data']) for r in db.execute('select * from markers')}
        genotypes={r['id']:json.loads(r['data']) for r in db.execute('select * from genotypes')}
        sets={r['id']:json.loads(r['data']) for r in db.execute('select * from genosets')}
    finally:db.close()
    from .assay import probe_warning_records
    probe_warnings=probe_warning_records(markers,genotypes)
    aliases=defaultdict(set);rejected_aliases=defaultdict(set)
    for marker,m in markers.items():
        if m['alias'] and marker in genome.by_id:aliases[m['alias']].add(marker)
        placement=m.get('public_placement')
        if placement and placement['build']==genome.build:
            for synonym in placement['synonyms']+[placement['canonical_id']]:
                if synonym and synonym!=marker and synonym in genome.by_id:
                    if all((r.chrom,r.pos)==(placement['chrom'],placement['pos']) for r in genome.by_id[synonym]):aliases[marker].add(synonym)
                    else:rejected_aliases[marker].add(synonym)
    findings=[];audit=[];calls={};evidence={};comp=str.maketrans('ACGT','TGCA')
    for marker in sorted(set(genome.by_id)&markers.keys() | aliases.keys() | rejected_aliases.keys()):
        m=markers.get(marker)
        if not m:continue
        ids=({marker} if marker in genome.by_id else set()) | aliases[marker]
        if not ids:
            audit.append(dict(id=marker,reason='documented_alias_coordinate_conflict',raw_ids=sorted(rejected_aliases[marker])));continue
        positions={(r.chrom,r.pos) for i in ids for r in genome.by_id.get(i,[])}
        reason=None;pos=m['positions'].get(genome.build)
        placement=m.get('public_placement');placement=placement if placement and placement['build']==genome.build else None
        verified_orientation=m.get('orientations_by_build',{}).get(genome.build)
        orientation=verified_orientation or m['orientation']
        own_snv={r.gt for i in ids for r in genome.by_id.get(i,[]) if r.gt!='--'}
        invariant=bool(placement and len(own_snv)==1 and set(next(iter(own_snv))) in ({'A','T'},{'C','G'}) and set(next(iter(own_snv)))==set(m.get('source_snv_alleles',[])) and set(next(iter(own_snv)))<=set(placement['alleles']))
        if invariant:orientation='plus'
        if m['alias']:continue  # Resolve chip probes through their explicit public rsID mapping.
        decoded=resolve_indel(genome,m,ids)
        if decoded and decoded['coding_status']=='provider_coding_unverified':
            literal=marker+'('+';'.join(sorted(decoded['raw_genotype'] if len(decoded['raw_genotype'])==2 else decoded['raw_genotype']*2)).lower()+')'
            if literal not in genotypes or positions!={(m['chrom'],pos)}:
                candidate_key=marker+'('+';'.join(decoded['source_alleles']).lower()+')'
                candidate=genotypes.get(candidate_key,{})
                audit.append(dict(id=marker,reason='provider_indel_coding_unverified',raw_ids=sorted(ids),
                    sequence_candidate=decoded,candidate_summary=candidate.get('summary',''),
                    explanation='The public sequence is known, but the provider-specific homozygous I/D allele mapping is not established. The tentative conversion is not a personal finding.'))
                continue
            # A literal source I/D page needs no sequence polarity assumption.
            decoded=None
        if not pos:reason='source_coordinate_unresolved'
        elif positions!={(m['chrom'],pos)}:reason='coordinate_or_alias_conflict'
        elif orientation not in {'plus','minus'}:reason='source_orientation_unresolved'
        elif (m['flip'+genome.build] or m['mapping_strand']=='-') and not verified_orientation and not invariant:reason='build_orientation_exception'
        if encoding_warning(m) and any(set(r.gt)<=set('ID') for i in ids for r in genome.by_id.get(i,[])):
            reason='source_indel_encoding_warning'
        if reason and not decoded:
            audit.append(dict(id=marker,reason=reason,raw_ids=sorted(ids)));continue
        status,gt,raw_ids=genome.locus(m['chrom'],pos) if pos else ('source_coordinate_unresolved','',sorted(ids))
        # A deletion assay and a base substitution can share a coordinate.
        # Only separate the assay types when IDs/explicit aliases identify this
        # variant; retain conflicts among all co-located base-substitution probes.
        own_rows=[r for i in ids for r in genome.by_id.get(i,[])]
        own_calls={r.gt for r in own_rows if r.gt!='--'}
        # A different rsID at the same base cannot supply this assay's no-call.
        # Only the exact ID or a documented public alias can provide evidence.
        if not own_calls:status='no_call'
        if status=='conflicting_probes' and own_calls and all(set(g)<=set('ACGT') for g in own_calls):
            snvs=[r for r in genome.by_position[(m['chrom'],pos)] if r.gt=='--' or set(r.gt)<=set('ACGT')]
            observed={r.gt for r in snvs if r.gt!='--'}
            if len(observed)==1 and not any(r.marker in genome.bad_ids for r in snvs):
                gt=next(iter(observed));status='called' if not m['chrom'].isdigit() or len(gt)==2 else 'unexpected_autosomal_ploidy'
        # Literal I/D source pages need no invented sequence conversion. Other
        # indels remain unresolved, even when a nearby genotype page exists.
        if status=='indel_code_unresolved' and genome.qc.get('provider')=='23andMe' and set(gt)<=set('ID'):
            key=marker+'('+';'.join(sorted(gt if len(gt)==2 else gt*2)).lower()+')'
            if key in genotypes and (not m['chrom'].isdigit() or len(gt)==2):status='called'
        if decoded:status,gt,raw_ids='called',decoded['raw_genotype'],sorted(ids)
        if status=='called' and placement and not set(gt)<=set(placement['alleles']):status='public_reference_allele_mismatch'
        definition=m.get('curated_snv_definition');target=definition['placements'].get(genome.build) if definition else None
        if status=='called' and target and not set(gt)<=set(target['ref']+target['alt']):status='public_reference_allele_mismatch'
        if status!='called':audit.append(dict(id=marker,reason=status,raw_ids=raw_ids));continue
        source_gt=gt.translate(comp) if orientation=='minus' else gt
        if len(source_gt)==1:source_gt*=2
        source_gt=tuple(sorted(source_gt));calls[marker]=source_gt
        if decoded:source_gt=decoded['source_alleles'];calls[marker]=source_gt
        flags=['published_array_quality_flag'] if (m['chrom'],pos) in flagged else []
        if m.get('coordinate_audit'):flags.append('source_coordinate_corrected_from_ClinVar')
        if placement:flags.append('public_SNV_placement_corroborated')
        if verified_orientation:flags.append('source_orientation_corroborated_from_public_alleles')
        if invariant:flags.append('strand_invariant_heterozygote')
        if rejected_aliases[marker]:flags.append('incompatible_alias_coordinates_excluded')
        if set(gt)<=set('ID'):flags.append('public_sequence_resolved_indel' if decoded else 'literal_source_indel_code_match')
        evidence[marker]=dict(raw_ids=sorted(ids),public_aliases=sorted(aliases[marker]),colocated_raw_ids=raw_ids,chrom=m['chrom'],pos=pos,raw_genotype=gt,source_genotype='('+';'.join(source_gt)+')',flags=flags,coordinate_audit=m.get('coordinate_audit'))
        cautions=[r for i in sorted(ids) if any(x.gt!='--' for x in genome.by_id.get(i,[])) for r in probe_warnings.get(i,[])]
        if cautions:evidence[marker]['probe_assay_warnings']=cautions;flags.append('input_probe_source_assay_warning')
        if placement:evidence[marker]['public_placement_audit']={**m['public_placement_audit'],**placement}
        if rejected_aliases[marker]:evidence[marker]['excluded_public_aliases']=sorted(rejected_aliases[marker])
        if target:
            flags.append('reviewed_public_sequence_definition')
            evidence[marker]['curated_sequence_evidence']=dict(**definition,build=genome.build,observed_target_alt_copies=gt.count(target['alt']),orientation=orientation)
        if decoded:
            evidence[marker]['sequence_resolution']=decoded
            if any((m['chrom'],p) in flagged for p in decoded['raw_positions']) and 'published_array_quality_flag' not in flags:flags.append('published_array_quality_flag')
        key=marker+'('+';'.join(source_gt).lower()+')';specific=genotypes.get(key)
        # A marker can have source metadata but no prose or genotype summary.
        findings.append(dict(id=marker,kind='variant',gene=m['gene'],summary=specific['summary'] if specific and specific['summary'] else 'Variant literature available; no genotype-specific summary in this snapshot.',magnitude=specific['magnitude'] if specific else 0,repute=specific['repute'] if specific else 'Unknown',genotype_text=specific['text'] if specific else '',variant_text=m['text'],pmids=sorted(set(m['pmids']+(specific['pmids'] if specific else []))),populations=m['populations'],source=specific['source'] if specific else m['source'],has_genotype_summary=bool(specific and specific['summary']),source_date=specific['retrieved'] if specific else m['retrieved'],**evidence[marker]))
        if target and gt.count(target['alt'])==0:
            findings[-1]['original_source_summary']=findings[-1]['summary']
            findings[-1]['summary']='Target alternate allele not observed at this site: '+definition['name']
    for row in audit:
        source=markers[row['id']]
        row.update(gene=source['gene'],source=source['source'],variant_context=source['text'])
    criteria={k:r['rule'] for k,r in sets.items() if r['rule']};outcomes,errors=evaluate_all(criteria,calls)
    partial=conditional_evidence(criteria,calls,outcomes)
    def deps(key,seen=None):
        seen=set() if seen is None else seen
        if key in seen:return set()
        seen.add(key);text=criteria.get(key,'');ids=set(re.findall(r'\b(?:rs|i)\d+',re.sub(r'(?m)^\s*#.*$','',text)))
        for child in re.findall(r'\bgs\d+',text):ids|=deps(child,seen)
        return ids
    genoset_audit=[]
    unresolved_by_id={a['id']:a['reason'] for a in audit}
    for key,r in sets.items():
        state=outcomes.get(key);dependencies=deps(key);missing=sorted(dependencies-calls.keys());flags=sorted({f for i in dependencies for f in evidence.get(i,{}).get('flags',[])})
        genoset_audit.append(dict(id=key,summary=r['summary'],source=r['source'],rule=r['rule'],status='matched' if state is True else 'not_matched' if state is False else 'unresolved',reason=errors.get(key,'rule_unavailable' if not r['rule'] else 'missing_or_unresolved_dependency' if state is None else ''),missing_markers=missing,missing_details=[dict(id=i,reason=unresolved_by_id.get(i,'not_assayed' if i not in genome.by_id else 'source_record_unavailable')) for i in missing],flags=flags))
        genoset_audit[-1]['blocking_markers']=unresolved_markers(key,criteria,calls,outcomes) if state is None else []
        if key in partial:genoset_audit[-1]['conditional_evidence']=partial[key]
        if state is True:findings.append(dict(id=key,kind='genoset',gene='',summary=r['summary'],magnitude=r['magnitude'],repute=r['repute'],genotype_text=r['text'],variant_text='',pmids=r['pmids'],populations=[],source=r['source'],source_date=r['rule_date'],rule_revision=r['revision'],dependencies=sorted(dependencies),dependency_calls={i:evidence[i] for i in dependencies if i in evidence},flags=flags,has_genotype_summary=bool(r['summary'])))
    findings.sort(key=lambda r:(-r['magnitude'],r['id']))
    return dict(findings=findings,unresolved=audit,genoset_audit=genoset_audit,counts=dict(findings=len(findings),variant_summaries=sum(r['kind']=='variant' and r['has_genotype_summary'] for r in findings),genosets=sum(r['kind']=='genoset' for r in findings),matched_loci=len(calls),unresolved_loci=len(audit)),reference=meta)
