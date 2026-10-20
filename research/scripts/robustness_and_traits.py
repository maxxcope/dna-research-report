"""Local sensitivity checks and a transparent historical six-SNP eye-colour model.

No network calls. Model coefficients are the original IrisPlex coefficients,
as reproduced in the archived UCSD teaching workbook and IU thesis Table 2.2.
Alleles are explicitly harmonized to the raw file's GRCh37 forward strand.
This is not the current Erasmus online calculator or a clinical prediction.
"""
import json
from pathlib import Path
from collections import Counter
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA

ROOT=Path(__file__).resolve().parents[1]
a=ROOT/'analysis'
d=np.load(a/'ancestry_genotypes.npz',allow_pickle=True) # own locally generated artifact
bim=pd.read_csv(a/'ancestry_reference.bim',sep=r'\s+',header=None)
results={}
for parity in (0,1):
    sel=bim[0].to_numpy()%2==parity
    for name,smask,k in [('global',np.ones(len(d['X']),bool),10),('europe',d['supers']=='EUR',6)]:
        x=d['X'][smask][:,sel];u=d['user'][sel];lab=d['labels'][smask]
        mu=np.nanmean(x,axis=0);p=mu/2;good=(p>.05)&(p<.95)
        sd=np.sqrt(2*p[good]*(1-p[good]))
        z=np.nan_to_num((x[:,good]-mu[good])/sd,nan=0)
        model=PCA(n_components=k,svd_solver='randomized',random_state=45,iterated_power=7)
        ref=model.fit_transform(z);up=model.transform(((u[good]-mu[good])/sd)[None,:])[0]
        near=np.argsort(np.linalg.norm(ref-up,axis=1))[:30]
        results[f'{name}_{"even" if parity==0 else "odd"}_chromosomes']={'markers':int(good.sum()),'nearest_30':dict(Counter(lab[near]))}

raw=json.loads((a/'raw_by_id.json').read_text())
parameters=[('rs12913832','A',-4.81,-1.79),('rs1800407','T',1.4,.87),('rs12896399','G',-.58,-.03),('rs16891982','C',-1.3,-.5),('rs1393350','A',.47,.27),('rs12203592','T',.7,.73)]
blue,other=3.94,.65
rows=[]
for rs,allele,b1,b2 in parameters:
    gt=raw[rs]['gt'];assert len(gt)==2 and set(gt)<=set('ACGT')
    dose=gt.count(allele);blue+=dose*b1;other+=dose*b2
    rows.append(dict(rsid=rs,raw_genotype=gt,effect_allele_forward=allele,dosage=dose,beta_blue=b1,beta_intermediate=b2))
weights=np.exp([blue,other,0]);probs=weights/weights.sum()
results['historical_irisplex']={'model':'Original six-SNP multinomial model, rounded published coefficients; not current web fit','inputs':rows,'intercepts':{'blue':3.94,'intermediate':.65,'brown_reference':0},'probabilities':dict(zip(['blue','intermediate','brown'],probs.tolist())),'caution':'Model estimates, not measured eye colour or individually calibrated confidence.'}
(a/'robustness_and_traits.json').write_text(json.dumps(results,indent=2,default=lambda x:int(x)))
print(json.dumps(results,indent=2,default=lambda x:int(x)))
