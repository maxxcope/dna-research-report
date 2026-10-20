"""Public reference preparation, deliberately without any DNA input argument."""
import argparse
from collections import defaultdict
import csv
import fcntl
from datetime import datetime, timezone
import json
import re
from pathlib import Path
import shutil
import signal
import subprocess
import tarfile
import urllib.request
import zipfile
from .downloads import download as transfer, atomic_json
from .references import digest, fetch_clinvar, import_gwas, metadata
from .knowledge import import_knowledge
from .scoring import score_rows

SNPEDIA='https://static.snpbrowser.com/snpedia.db'
VCF='https://bochet.gcc.biostat.washington.edu/beagle/1000_Genomes_phase3_v5a/b37.vcf/chr{ch}.1kg.phase3.v5a.vcf.gz'

def download(url,path,manifest):
    path=Path(path)
    previous=manifest.get(str(path.name))
    if previous and previous.get('url')!=url:
        raise ValueError('Cached reference source changed: '+path.name+'. Prepare a new cache directory.')
    transfer(url,path,expected=previous['sha256'] if previous else None)
    manifest[path.name]=dict(url=url,sha256=digest(path),bytes=path.stat().st_size)
    return path

def extract_files(archive,destination,names):
    """Extract named regular files only; never archive paths or symlinks."""
    destination=Path(destination);destination.mkdir(parents=True,exist_ok=True)
    if tarfile.is_tarfile(archive):
        with tarfile.open(archive) as source:
            for name in names:
                entry=source.getmember(name)
                if not entry.isfile():raise ValueError('Archive member is not a regular file')
                with source.extractfile(entry) as inp,(destination/Path(name).name).open('wb') as out:shutil.copyfileobj(inp,out)
    else:
        with zipfile.ZipFile(archive) as source:
            for name in names:
                with source.open(name) as inp,(destination/Path(name).name).open('wb') as out:shutil.copyfileobj(inp,out)

def planned_windows(score):
    groups=defaultdict(list)
    for r in score:groups[r['hm_chr']].append((max(1,int(r['hm_pos'])-1000000),int(r['hm_pos'])+1000000))
    out=[]
    for ch,intervals in sorted(groups.items(),key=lambda x:int(x[0])):
        merged=[]
        for start,end in sorted(intervals):
            if merged and start<=merged[-1][1]:merged[-1][1]=max(end,merged[-1][1])
            else:merged.append([start,end])
        out.extend(dict(name=f'chr{ch}_{a}_{b}',chrom=ch,start=a,end=b) for a,b in merged)
    return out

def prepare_windows(windows,root,samples,manifest,checkpoint=lambda:None):
    import pysam
    groups=defaultdict(list)
    for w in windows:groups[w['chrom']].append(w)
    out=[]
    for ch,group in groups.items():
        url=VCF.format(ch=ch);remote=None
        try:
            for w in group:
                p=root/'windows'/(w['name']+'.vcf.gz');p.parent.mkdir(exist_ok=True)
                if not p.exists():
                    if remote is None:
                        # HTSlib otherwise saves remote indexes in the working
                        # directory. Keep every downloaded reference in the cache.
                        index=download(url+'.tbi',root/'remote-indexes'/f'remote-chr{ch}.vcf.gz.tbi',manifest)
                        checkpoint()
                        remote=pysam.VariantFile(url,index_filename=str(index));remote.subset_samples(samples)
                    print('Preparing public window '+w['name'],flush=True);part=Path(str(p)+'.part')
                    try:
                        with pysam.VariantFile(str(part),'wz',header=remote.header) as dst:
                            for rec in remote.fetch(ch,w['start']-1,w['end']):
                                if len(rec.alleles)>1 and all(len(a)==1 and a in 'ACGT' for a in rec.alleles):dst.write(rec)
                        part.replace(p)
                    finally:part.unlink(missing_ok=True)
                checksum=digest(p);previous=manifest.get(w['name'])
                if previous and previous['sha256']!=checksum:raise ValueError('Cached window checksum changed')
                if not Path(str(p)+'.tbi').exists():pysam.tabix_index(str(p),preset='vcf',force=True)
                with pysam.VariantFile(p) as vf:
                    if list(vf.header.samples)!=samples:raise ValueError('Reference window has unexpected sample order')
                    if any(str(rec.chrom)!=ch or not w['start']<=rec.pos<=w['end'] for rec in vf):raise ValueError('Reference window has unexpected coordinates')
                manifest[w['name']]=dict(url=url,sha256=checksum,start=w['start'],end=w['end'],reference_samples=len(samples))
                checkpoint()
                out.append({**w,'file':str(p),'map':str(root/'maps'/f'plink.chr{ch}.GRCh37.map'),'source':url,'sha256':checksum})
                print(f"Public windows: {len(out)}/{len(windows)} checked",flush=True)
        finally:
            if remote:remote.close()
    return out

def prepare(root,plink,java,progress=print):
    root=Path(root).expanduser().resolve();root.mkdir(parents=True,exist_ok=True)
    # Different application installations can share the same reference cache.
    with (root/'.prepare.lock').open('a+') as lock:
        try:fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        except BlockingIOError:raise ValueError('Reference preparation is already running for this folder.') from None
        try:return _prepare(root,plink,java,progress)
        finally:fcntl.flock(lock,fcntl.LOCK_UN)

def _prepare(root,plink,java,progress=print):
    root=Path(root).expanduser().resolve();root.mkdir(parents=True,exist_ok=True)
    plink=str(Path(shutil.which(str(plink)) or plink).resolve());java=str(Path(shutil.which(str(java)) or java).resolve())
    subprocess.run([plink,'--version'],check=True,capture_output=True)
    subprocess.run([java,'-version'],check=True,capture_output=True)
    manifest_path=root/'downloads.json';manifest=json.loads(manifest_path.read_text()) if manifest_path.exists() else {}
    if not isinstance(manifest,dict) or any(not isinstance(r,dict) or not re.fullmatch(r'[0-9a-f]{64}',str(r.get('sha256',''))) for r in manifest.values()):
        raise ValueError('The download inventory is malformed. Preserve this cache and choose a new reference folder.')
    def get(url,name):
        path=download(url,root/name,manifest);atomic_json(manifest_path,manifest);return path
    chain=get('https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz','hg19ToHg38.over.chain.gz')
    source=get(SNPEDIA,'snpedia.db');criteria=root/'genoset-criteria.json'
    from .fetch_rules import download as fetch_rules
    fetch_rules(source,criteria)
    prepared=metadata(root/'snpedia.sqlite') if (root/'snpedia.sqlite').exists() else {}
    if not (root/'clinvar37.sqlite').exists():
        from .references import import_clinvar
        url='https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/clinvar.vcf.gz'
        progress('Preparing ClinVar index')
        import_clinvar(get(url,'clinvar37.vcf.gz'),root/'clinvar37.sqlite','37',url,progress)
    from .source_enrichment import download as enrich_source
    placements=root/'ensembl37.json';enrich_source(source,criteria,chain,placements,progress)
    from .curation import definitions
    if prepared.get('importer_version')!='6' or prepared.get('curated_metadata_sha256')!=definitions()[1] or prepared.get('source_sha256')!=digest(source) or prepared.get('criteria_sha256')!=digest(criteria) or prepared.get('chain_sha256')!=digest(chain) or prepared.get('coordinate_reference_sha256')!=digest(root/'clinvar37.sqlite') or prepared.get('placement_reference_sha256')!=digest(placements):
        progress('Building SNPedia interpretation index')
        import_knowledge(source,criteria,chain,root/'snpedia.sqlite',root/'clinvar37.sqlite',placements,progress)
        progress('SNPedia interpretation index complete')
    gw=get('https://ftp.ebi.ac.uk/pub/databases/gwas/releases/latest/gwas-catalog-associations_ontology-annotated-full.zip','gwas.zip')
    if not (root/'gwas.sqlite').exists():import_gwas(gw,root/'gwas.sqlite',manifest[gw.name]['url'],progress)
    clusters=get('https://zenodo.org/records/5047472/files/the_list.tsv.gz','chip_clusters.tsv.gz')
    bad=get('https://zenodo.org/records/5047472/files/badalleles.tsv.gz','badalleles.tsv.gz')
    panel=get('https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel','1000g_panel.tsv')
    tcag=get('https://tcag.ca/documents/tools/omni25_indep.tar.gz','omni25_indep.tar.gz')
    if not all((root/'omni25'/('indep.'+s)).exists() for s in ['bed','bim','fam']):extract_files(tcag,root/'omni25',['indep.bed','indep.bim','indep.fam'])
    for name in ['hgdp650.qc.hg19.bed','hgdp650.qc.hg19.bim','hgdp650.qc.hg19.fam','hgdp650_id_pop_coords.txt']:get('https://zenodo.org/records/17375732/files/'+name,name)
    tree=get('https://www.familytreedna.com/public/y-dna-haplotree/get','ftdna_tree.json')
    beagle=get('https://faculty.washington.edu/browning/beagle/beagle.27Feb25.75f.jar','beagle.jar')
    maps=get('https://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/plink.GRCh37.map.zip','maps.zip')
    mapnames=[f'plink.chr{i}.GRCh37.map' for i in range(1,23)]
    if not all((root/'maps'/n).exists() for n in mapnames):extract_files(maps,root/'maps',mapnames)
    hg=get('https://github.com/genepi/haplogrep3/releases/download/v3.3.2/haplogrep3-3.3.2-linux.tar.gz','haplogrep3.tar.gz');hgdir=root/'haplogrep3'
    if not (hgdir/'haplogrep3.jar').exists():extract_files(hg,hgdir,['haplogrep3.jar','haplogrep3.yaml'])
    if not (hgdir/'trees/phylotree-rcrs/17.3/rcrs.fasta').exists():
        subprocess.run([java,'-jar',str(hgdir/'haplogrep3.jar'),'install-tree','phylotree-rcrs@17.3'],cwd=hgdir,check=True)
    scores=[]
    for sid,name in [('PGS000034','BMI'),('PGS000019','Coronary artery disease'),('PGS000297','Height')]:
        path=get(f'https://ftp.ebi.ac.uk/pub/databases/spot/pgs/scores/{sid}/ScoringFiles/Harmonized/{sid}_hmPOS_GRCh37.txt.gz',sid+'.txt.gz')
        scores.append(dict(id=sid,name=name,path=str(path),build='37',source=f'https://www.pgscatalog.org/score/{sid}/'))
    samples=[r['sample'] for r in csv.DictReader(panel.open(),delimiter='\t') if r['super_pop']=='EUR']
    windows=planned_windows(score_rows(scores[0]['path']))
    checkpoint=lambda:atomic_json(manifest_path,manifest)
    windows=prepare_windows(windows,root,samples,manifest,checkpoint)
    phase_window=prepare_windows([dict(name='chr6_17100000_19100000',chrom='6',start=17100000,end=19100000)],root,samples,manifest,checkpoint)[0]
    profile=dict(schema=1,knowledge=str(root/'snpedia.sqlite'),clinvar=str(root/'clinvar37.sqlite'),gwas=str(root/'gwas.sqlite'),quality=dict(clusters=str(clusters),badalleles=str(bad)),ancestry=[dict(name='1000 Genomes',prefix=str(root/'omni25/indep'),labels=str(panel),plink=plink,build='37',source='https://tcag.ca/tools/1000genomes.html'),dict(name='HGDP',prefix=str(root/'hgdp650.qc.hg19'),labels=str(root/'hgdp650_id_pop_coords.txt'),plink=plink,build='37',sample_column='id',population_column='population',region_column='region',source='https://zenodo.org/records/17375732')],frequencies=dict(prefix=str(root/'omni25/indep'),labels=str(panel),source='https://tcag.ca/tools/1000genomes.html'),paternal=dict(tree=str(tree),chain=str(chain),source='https://www.familytreedna.com/public/y-dna-haplotree/get'),maternal=dict(directory=str(hgdir),java=java,tree='phylotree-rcrs@17.3'),scores=scores,imputation=dict(id='PGS000034',name='BMI research score',score=scores[0]['path'],windows=windows,java=java,beagle=str(beagle),build='37',workers=2,validation_windows=2,source=scores[0]['source']),phase=dict(name='TPMT target phase',java=java,beagle=str(beagle),build='37',window=phase_window,sites=[dict(id='rs1800460',chrom='6',pos=18139228,effect_allele='T'),dict(id='rs1142345',chrom='6',pos=18130918,effect_allele='C')],source='https://cpicpgx.org/guidelines/guideline-for-thiopurines-and-tpmt/'))
    # Relative data paths make the cache movable; installed executables stay explicit.
    def relative(value):
        if isinstance(value,dict):return {k:relative(v) for k,v in value.items()}
        if isinstance(value,list):return [relative(v) for v in value]
        if isinstance(value,str) and value.startswith(str(root)+'/'):return str(Path(value).relative_to(root))
        return value
    manifest_path.write_text(json.dumps(manifest,indent=2)+'\n')
    used=set()
    def collect(value,key=''):
        if isinstance(value,dict):
            for k,v in value.items():collect(v,k)
        elif isinstance(value,list):
            for v in value:collect(v,key)
        elif isinstance(value,str):
            if key=='prefix':used.update(Path(value+'.'+ext) for ext in ['bed','bim','fam'])
            elif key=='directory':used.update(p for p in Path(value).rglob('*') if p.is_file())
            elif key in {'clusters','badalleles','labels','plink','tree','chain','java','beagle','score','file','map','path','knowledge','clinvar','gwas'} and Path(value).is_file():used.add(Path(value))
    collect(profile)
    used.add(placements)
    profile['reference_files']=[dict(file=str(p),sha256=digest(p)) for p in sorted(used)]
    profile['reference_prepared_utc']=datetime.now(timezone.utc).isoformat()
    tmp=root/'pipeline.json.tmp';tmp.write_text(json.dumps(relative(profile),indent=2)+'\n');tmp.replace(root/'pipeline.json')
    progress('Full public reference profile ready: '+str(root/'pipeline.json'))

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--refs',type=Path,required=True);parser.add_argument('--plink',default='plink');parser.add_argument('--java',default='java')
    args=parser.parse_args()
    def cancelled(signum,frame):raise KeyboardInterrupt()
    for sig in (signal.SIGTERM,signal.SIGHUP):signal.signal(sig,cancelled)
    try:prepare(args.refs,args.plink,args.java)
    except KeyboardInterrupt:
        parser.exit(130,'Reference preparation cancelled; completed files are kept.\n')
    except (OSError,ValueError) as error:
        parser.exit(1,'Reference preparation stopped: '+str(error)+'\n')

if __name__=='__main__':main()
