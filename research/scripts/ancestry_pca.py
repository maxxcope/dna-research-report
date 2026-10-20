"""Reference-only PCA and projection; allele count convention verified against PLINK."""
import os
os.environ.setdefault('MPLCONFIGDIR','/tmp/genetics-methods-mpl')
from pathlib import Path
import numpy as np,pandas as pd,json
from collections import Counter
from sklearn.decomposition import PCA
from sklearn.neighbors import KNeighborsClassifier
from sklearn.model_selection import train_test_split
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from inventory import ROOT

prefix=ROOT/'analysis/ancestry_reference'
bim=pd.read_csv(str(prefix)+'.bim',sep=r'\s+',header=None,names=['ch','id','cm','pos','a1','a2'])
fam=pd.read_csv(str(prefix)+'.fam',sep=r'\s+',header=None,names=['fid','iid','pat','mat','sex','pheno'])
panel=pd.read_csv(ROOT/'references/1000g_panel.tsv',sep='\t').set_index('sample')
labels=panel.loc[fam.iid,'pop'].to_numpy();supers=panel.loc[fam.iid,'super_pop'].to_numpy()
n=len(fam);m=len(bim);stride=(n+3)//4
with open(str(prefix)+'.bed','rb') as f:
    assert f.read(3)==bytes([0x6c,0x1b,0x01]);buf=np.frombuffer(f.read(),dtype=np.uint8).reshape(m,stride)
codes=((buf[:,:,None] >> np.array([0,2,4,6],dtype=np.uint8)) & 3).reshape(m,-1)[:,:n]
# PLINK BED 00=A1/A1, 01=missing, 10=A1/A2, 11=A2/A2.
X=np.array([2,np.nan,1,0],dtype=np.float32)[codes].T.copy()
mapping={r['id']:r for r in json.loads((ROOT/'analysis/ancestry_marker_mapping.json').read_text())}
user=np.array([mapping[r.id]['gt'].count(r.a1) for r in bim.itertuples()],dtype=np.float32)
np.savez_compressed(ROOT/'analysis/ancestry_genotypes.npz',X=X,user=user,labels=labels,supers=supers)

def run(mask,name,k):
    xx=X[mask]; labs=labels[mask]; sups=supers[mask]
    mu=np.nanmean(xx,axis=0);af=mu/2
    good=(af>.05)&(af<.95)
    sd=np.sqrt(2*af[good]*(1-af[good]))
    z=np.nan_to_num((xx[:,good]-mu[good])/sd,nan=0)
    u=((user[good]-mu[good])/sd).reshape(1,-1)
    pca=PCA(n_components=k,svd_solver='randomized',random_state=42,iterated_power=7)
    pcs=pca.fit_transform(z);up=pca.transform(u)[0]
    distances=np.linalg.norm(pcs[:,:k]-up[:k],axis=1)
    near=np.argsort(distances)[:30]
    centroids={l:float(np.linalg.norm(pcs[labs==l,:min(k,4)].mean(axis=0)-up[:min(k,4)])) for l in set(labs)}
    # Hold out reference samples, then fit PCA exclusively to the training subset.
    tr,te=train_test_split(np.arange(len(xx)),test_size=.2,random_state=42,stratify=labs)
    muv=np.nanmean(xx[tr][:,good],axis=0);afv=muv/2
    sdv=np.sqrt(np.maximum(2*afv*(1-afv),.01))
    zz=np.nan_to_num((xx[:,good]-muv)/sdv,nan=0)
    pv=PCA(n_components=k,svd_solver='randomized',random_state=43,iterated_power=7)
    a=pv.fit_transform(zz[tr]);b=pv.transform(zz[te])
    knn=KNeighborsClassifier(n_neighbors=15).fit(a,labs[tr]);pred=knn.predict(b)
    acc=float((pred==labs[te]).mean())
    super_lookup=dict(zip(labs,sups))
    super_acc=float(np.mean([super_lookup[p]==s for p,s in zip(pred,sups[te])]))
    nearest=Counter(labs[near]);print(name,'markers',good.sum(),'nearest30',nearest,'CV pop',acc,'CV super',super_acc,flush=True)
    out=dict(markers=int(good.sum()),references=len(xx),components=k,explained_variance_ratio=pca.explained_variance_ratio_.tolist(),user_pc=up.tolist(),
        nearest_30_population_counts={str(l):int(c) for l,c in nearest.items()},centroid_distances_first4=centroids,
        heldout_population_accuracy=acc,heldout_superpopulation_accuracy=super_acc)
    df=pd.DataFrame(pcs,columns=[f'PC{i+1}' for i in range(k)]);df.insert(0,'population',labs);df.insert(0,'sample',fam.iid.to_numpy()[mask])
    df.to_csv(ROOT/f'analysis/{name}_reference_pcs.tsv',sep='\t',index=False)
    for axes in [(0,1),(2,3)]:
        fig,ax=plt.subplots(figsize=(9,6.5),layout='constrained')
        ll=sups if name=='global' else labs
        colors=plt.colormaps['tab10']
        for i,l in enumerate(sorted(set(ll))):
            w=ll==l;ax.scatter(pcs[w,axes[0]],pcs[w,axes[1]],s=13,alpha=.55,label=l,color=colors(i))
        ax.scatter(up[axes[0]],up[axes[1]],s=200,marker='*',edgecolor='black',color='#ef3454',label='Your projected sample',zorder=5)
        ax.set(xlabel=f'PC{axes[0]+1}',ylabel=f'PC{axes[1]+1}',title=f'{"Global" if name=="global" else "European reference"} genetic similarity • {int(good.sum()):,} SNPs')
        ax.legend(fontsize=9,frameon=False);ax.spines[['top','right']].set_visible(False)
        fig.savefig(ROOT/f'output/{name}_pca_{axes[0]+1}_{axes[1]+1}.png',dpi=160);plt.close(fig)
    return out

summary=dict(global_pca=run(np.ones(n,dtype=bool),'global',10),europe_pca=run(supers=='EUR','europe',6))
(ROOT/'analysis/ancestry_summary.json').write_text(json.dumps(summary,indent=2))
