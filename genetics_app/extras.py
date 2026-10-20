"""Input-driven research analyses. Public references are read, never modified."""
from collections import Counter, defaultdict
import csv
import gzip
import json
from pathlib import Path
import subprocess
import os
import signal
import threading
from .knowledge import Chain


_commands = set()
_command_lock = threading.Lock()
_cancelled = False


def reset_commands():
    global _cancelled
    with _command_lock:
        if _commands:
            raise ValueError('An analysis is already running in this process.')
        _cancelled = False


def _stop_command(process):
    if process.poll() is not None:
        return
    try:
        if os.name == 'posix':
            os.killpg(process.pid, signal.SIGTERM)
        else:
            process.terminate()
        try:
            process.wait(timeout=2)
        except subprocess.TimeoutExpired:
            if os.name == 'posix':
                os.killpg(process.pid, signal.SIGKILL)
            else:
                process.kill()
            process.wait()
    except ProcessLookupError:
        pass


def cancel_commands():
    global _cancelled
    with _command_lock:
        _cancelled = True
        active = list(_commands)
    for process in active:
        _stop_command(process)


def run_command(command, cwd, timeout=900):
    with _command_lock:
        if _cancelled:
            raise InterruptedError('Analysis cancelled.')
        process = subprocess.Popen([str(x) for x in command], cwd=cwd,
                                   stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                   text=True, start_new_session=os.name == 'posix')
        _commands.add(process)
    try:
        stdout, stderr = process.communicate(timeout=timeout)
        if process.returncode:
            raise ValueError('Research tool failed: ' + Path(str(command[0])).name + '; ' + stderr[-600:])
        return stdout
    except BaseException:
        _stop_command(process)
        raise
    finally:
        for stream in (process.stdout, process.stderr):
            if stream is not None:
                stream.close()
        with _command_lock:
            _commands.discard(process)


def quality(genome,config):
    if genome.build!='37':return dict(status='unsupported_build',explanation='This GenomePrep reference uses GRCh37.'),set()
    bp=set(genome.by_position);totals=Counter();matches=Counter()
    with gzip.open(config['clusters'],'rt') as stream:
        for line in stream:
            site,groups=line.strip().split('\t');ch,pos=site.split(':');hit=(ch,int(pos)) in bp
            for group in groups.split(','):totals[group]+=1;matches[group]+=hit
    ranked=sorted([dict(cluster=k,overlap=matches[k],reference_sites=n,jaccard=matches[k]/(len(bp)+n-matches[k])) for k,n in totals.items()],key=lambda r:-r['jaccard'])
    if not ranked or ranked[0]['jaccard']<.5:return dict(status='insufficient_coverage',explanation='No sufficiently similar array-content cluster; no chip quality flags applied.',clusters=ranked),set()
    best=ranked[0]['cluster'];flagged=set()
    with gzip.open(config['badalleles'],'rt') as stream:
        for line in stream:
            cluster,sites=line.strip().split('\t')
            if cluster==best:
                for site in sites.split(','):
                    ch,pos=site.split(':');key=(ch,int(pos))
                    if key in bp:flagged.add(key)
    return dict(status='completed',selected_cluster=best,clusters=ranked,flagged_positions=len(flagged),explanation=f'{len(flagged):,} positions overlap published quality flags for the best-matching array cluster. These flags qualify findings; they do not prove any individual call wrong.',source='https://pmc.ncbi.nlm.nih.gov/articles/PMC8267563/'),flagged


def bed_matrix(prefix,keep_ids=None):
    import numpy as np
    prefix=str(prefix);fam=[l.split() for l in Path(prefix+'.fam').read_text().splitlines()];n=len(fam);stride=(n+3)//4
    records=[l.split() for l in Path(prefix+'.bim').read_text().splitlines()]
    with open(prefix+'.bed','rb') as f:
        if f.read(3)!=bytes([108,27,1]):raise ValueError('Expected SNP-major PLINK BED.')
    matrix=np.memmap(prefix+'.bed',dtype=np.uint8,mode='r',offset=3,shape=(len(records),stride))
    idx=[i for i,r in enumerate(records) if keep_ids is None or r[1] in keep_ids]
    if not idx:return np.empty((n,0),dtype=np.float32),[],fam
    buf=np.asarray(matrix[idx]);codes=((buf[:,:,None]>>np.array([0,2,4,6],dtype=np.uint8))&3).reshape(len(idx),-1)[:,:n]
    return np.array([2,np.nan,1,0],dtype=np.float32)[codes].T, [records[i] for i in idx],fam


def ancestry(genome,config,work,flagged):
    import numpy as np
    from sklearn.decomposition import PCA
    from sklearn.model_selection import train_test_split
    from sklearn.neighbors import KNeighborsClassifier
    from threadpoolctl import threadpool_limits
    if genome.build!=config.get('build','37'):return dict(status='unsupported_build',explanation='Ancestry panel assembly differs from the raw export.')
    prefix=config['prefix'];mapping={};ids=[]
    for line in open(str(prefix)+'.bim'):
        ch,marker,cm,pos,a1,a2=line.split();pos=int(pos)
        if ch not in {str(i) for i in range(1,23)} or (ch,pos) in flagged:continue
        if len(a1)!=1 or len(a2)!=1 or a1==a2 or set(a1+a2) in (set('AT'),set('CG')):continue
        state,gt,_=genome.locus(ch,pos)
        if state!='called' or len(gt)!=2 or not set(gt)<=set(a1+a2):continue
        mapping[marker]=gt;ids.append(marker)
    if len(ids)<500:return dict(status='insufficient_coverage',markers=len(ids),explanation='Fewer than 500 compatible autosomal reference markers.')
    work.mkdir(parents=True,exist_ok=True);extract=work/'extract.txt';extract.write_text('\n'.join(ids)+'\n')
    run_command([config['plink'],'--bfile',prefix,'--extract',extract,'--maf','.05','--geno','.02','--indep-pairwise','200','20','.2','--out',work/'pruned'],work)
    pruned=(work/'pruned.prune.in').read_text().splitlines()
    # Deterministic thinning bounds runtime; PLINK has already pruned LD.
    if len(pruned)>15000:pruned=pruned[::(len(pruned)+14999)//15000]
    X,bim,fam=bed_matrix(prefix,set(pruned))
    labels=list(csv.DictReader(open(config['labels']),delimiter='\t'))
    lookup={r[config.get('sample_column','sample')]:r for r in labels};idx=[i for i,r in enumerate(fam) if r[1] in lookup]
    X=X[idx];labs=np.array([lookup[fam[i][1]][config.get('population_column','pop')] for i in idx]);regions=np.array([lookup[fam[i][1]][config.get('region_column','super_pop')] for i in idx])
    if len(idx)<30 or X.shape[1]<100:return dict(status='insufficient_coverage',explanation='Too few labelled reference samples or pruned markers.')
    user=np.array([mapping[r[1]].count(r[4]) for r in bim],dtype=float)
    with threadpool_limits(limits=2):
        mu=np.nanmean(X,axis=0);af=mu/2;keep=(af>.05)&(af<.95);sd=np.sqrt(2*af[keep]*(1-af[keep]));Z=np.nan_to_num((X[:,keep]-mu[keep])/sd);u=(user[keep]-mu[keep])/sd
        k=min(10,len(idx)-2,Z.shape[1]-1);pca=PCA(n_components=k,svd_solver='randomized',random_state=2026);pcs=pca.fit_transform(Z);up=pca.transform(u[None])[0]
        distances=np.linalg.norm(pcs-up,axis=1);near=np.argsort(distances)[:30]
        centers=sorted([dict(population=str(l),distance=float(np.linalg.norm(pcs[labs==l].mean(axis=0)-up)),samples=int(sum(labs==l))) for l in set(labs)],key=lambda r:r['distance'])
        # Holdout PCA is fit only on training samples, with training-only scaling.
        counts=Counter(labs);eligible=np.array([i for i,l in enumerate(labs) if counts[l]>=5]);validation=None
        if len(eligible)>50:
            tr,te=train_test_split(eligible,test_size=.25,stratify=labs[eligible],random_state=2026)
            means=np.nanmean(X[tr],axis=0);freq=means/2;ok=(freq>.05)&(freq<.95);scale=np.sqrt(2*freq[ok]*(1-freq[ok]));zz=np.nan_to_num((X[:,ok]-means[ok])/scale)
            model=PCA(n_components=min(k,len(tr)-2,zz.shape[1]-1),svd_solver='randomized',random_state=2026);train=model.fit_transform(zz[tr]);test=model.transform(zz[te]);pred=KNeighborsClassifier(n_neighbors=min(7,len(tr))).fit(train,labs[tr]).predict(test)
            validation=dict(test_samples=len(te),population_accuracy=float(np.mean(pred==labs[te])))
    return dict(status='completed',name=config['name'],markers=int(keep.sum()),reference_samples=len(idx),nearest_30_population_counts=dict(Counter(map(str,labs[near]))),centroids=centers[:12],projection=up.tolist(),points=[dict(x=float(p[0]),y=float(p[1]),population=str(l),region=str(r)) for p,l,r in zip(pcs,labs,regions)],variance=pca.explained_variance_ratio_.tolist(),validation=validation,explanation='Genetic similarity to sampled reference populations, not ancestry percentages or nationality. The reference-only holdout accuracy measures population-label separation, not confidence in your ancestry.',source=config['source'])


def frequencies(genome,config,findings):
    import numpy as np
    selected={r['id'] for r in findings if r['kind']=='variant' and r['has_genotype_summary']}
    X,bim,fam=bed_matrix(config['prefix'],selected)
    labels={r['sample']:r for r in csv.DictReader(open(config['labels']),delimiter='\t')};out=[]
    for j,row in enumerate(bim):
        ch,marker,cm,pos,a1,a2=row;state,gt,_=genome.locus(ch,int(pos))
        if state!='called' or len(gt)!=2 or len(a1)!=1 or len(a2)!=1 or not set(gt)<=set(a1+a2):continue
        for region in sorted({v['super_pop'] for v in labels.values()}):
            idx=[i for i,r in enumerate(fam) if r[1] in labels and labels[r[1]]['super_pop']==region];v=X[idx,j];v=v[np.isfinite(v)]
            if len(v):out.append(dict(id=marker,population=region,matching_genotype=int(sum(v==gt.count(a1))),called_reference_samples=len(v),fraction=float(np.mean(v==gt.count(a1)))))
    return dict(status='completed',rows=out,explanation='Measured genotype proportions in the supplied public reference samples; not exact worldwide prevalence.',source=config['source'])


def y_lineage(genome,config):
    tree=json.loads(Path(config['tree']).read_text());nodes={int(k):v for k,v in tree['allNodes'].items()}
    for n in tree['roots']:nodes[n['haplogroupId']]=n
    chain=Chain(config['chain']) if genome.build=='37' else None
    bases={};excluded=0
    for (ch,pos),rows in genome.by_position.items():
        if ch!='Y':continue
        state,gt,_=genome.locus(ch,pos)
        if state!='called' or len(set(gt))!=1:excluded+=1;continue
        mapped=chain.map('Y',pos) if chain else (pos,'+')
        if not mapped:excluded+=1;continue
        base=gt[0].translate(str.maketrans('ACGT','TGCA')) if mapped[1]=='-' else gt[0]
        bases.setdefault(mapped[0],set()).add(base)
    if len(bases)<10:return dict(status='insufficient_coverage',explanation='Too few callable Y-chromosome positions; no paternal lineage assigned.')
    hits=defaultdict(list)
    for nid,n in nodes.items():
        for v in n.get('variants',[]):
            called=bases.get(v.get('position'))
            if not called or len(called)!=1 or len(v.get('ancestral',''))!=1 or len(v.get('derived',''))!=1:continue
            base=next(iter(called));state='derived' if base==v['derived'] else 'ancestral' if base==v['ancestral'] else 'other'
            hits[nid].append(dict(marker=v['variant'],position=v['position'],state=state))
    def path(nid):
        pp=[];seen=set()
        while nid in nodes and nid not in seen:seen.add(nid);pp.append(nid);nid=nodes[nid].get('parentId')
        return pp[::-1]
    candidates=[]
    for nid,rr in list(hits.items()):
        if not any(r['state']=='derived' for r in rr):continue
        pp=path(nid);allhits=[r for p in pp for r in hits.get(p,[])];der=sum(r['state']=='derived' for r in allhits);anc=sum(r['state']=='ancestral' for r in allhits)
        candidates.append(dict(lineage=nodes[nid]['name'],derived=der,ancestral_conflicts=anc,score=der-5*anc,terminal_support=sum(r['state']=='derived' for r in rr),path=[dict(lineage=nodes[p]['name'],calls=hits[p]) for p in pp if hits[p]]))
    candidates.sort(key=lambda x:(-x['score'],-x['derived'],-x['terminal_support'],x['lineage']))
    return dict(status='completed',tree_date=tree.get('publishedDate'),called_positions=len(bases),excluded_positions=excluded,candidates=candidates[:5],explanation='Candidate paternal branches ranked by observed derived support minus five times ancestral conflicts along the path. Unassayed descendants are not assigned; single-marker terminal branches remain tentative. This is a transparent tree comparison, not a calibrated probability.',source=config['source'])


def mito(genome,config,work):
    base=Path(config['directory']);tree=config.get('tree','phylotree-rcrs@17.3');tdir=base/'trees'/tree.replace('@','/');ref=''.join(l.strip() for l in (tdir/'rcrs.fasta').read_text().splitlines() if not l.startswith('>'))
    called={}
    for (ch,pos),rows in genome.by_position.items():
        if ch!='MT':continue
        state,gt,_=genome.locus(ch,pos)
        if state=='called' and len(set(gt))==1 and 1<=pos<=len(ref):called[pos]=gt[0]
    if len(called)<10:return dict(status='insufficient_coverage',explanation='Too few callable mitochondrial positions; no maternal lineage assigned.')
    work.mkdir(parents=True,exist_ok=True);hsd=work/'sample.hsd';hsd.write_text('LOCAL_SAMPLE\t"'+';'.join(map(str,sorted(called)))+'"\t\t'+'\t'.join(f'{p}{a}' for p,a in sorted(called.items()) if a!=ref[p-1])+'\n')
    output=work/'haplogrep.txt'
    run_command([config['java'],'-jar',base/'haplogrep3.jar','classify','--input',hsd,'--tree',tree,'--output',output,'--hits','5','--extend-report'],base)
    rows=list(csv.DictReader(output.open(),delimiter='\t'))
    return dict(status='completed',tree=tree,called_positions=len(called),results=rows,explanation='HaploGrep compares only the assayed mitochondrial positions with the chosen tree. Sparse array coverage can leave several maternal branches compatible.',source='https://haplogrep.readthedocs.io/')
