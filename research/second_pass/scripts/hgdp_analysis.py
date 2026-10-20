"""Broader reference ancestry comparison, with projection/joint-PCA sensitivity."""
import os
os.environ.setdefault('MPLCONFIGDIR','/tmp/genetics-methods-mpl')
from pathlib import Path
import json,collections,subprocess,sys
import numpy as np,pandas as pd
from sklearn.decomposition import PCA
from sklearn.model_selection import train_test_split
from sklearn.neighbors import KNeighborsClassifier
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[2];A=ROOT/'second_pass/analysis';R=ROOT/'second_pass/references';O=ROOT/'second_pass/output';L=ROOT/'second_pass/logs'
PLINK=ROOT/'references/plink/plink';pref=R/'hgdp650.qc.hg19'
cols=['ch','id','cm','pos','a1','a2']
def prepare():
 raw=json.loads((ROOT/'analysis/raw_by_id.json').read_text());bp=collections.defaultdict(list)
 for r in raw.values():
  if r['chrom'].isdigit() and 1<=int(r['chrom'])<=22:bp[(int(r['chrom']),r['pos'])].append(r)
 bim=pd.read_csv(str(pref)+'.bim',sep=r'\s+',header=None,names=cols);mapping=[];counts=collections.Counter()
 for r in bim.itertuples():
  if (r.ch,r.pos) not in bp:continue
  counts['overlapping_positions']+=1
  if len(str(r.a1))!=1 or len(str(r.a2))!=1 or not set(r.a1+r.a2)<=set('ACGT'):continue
  if set(r.a1+r.a2) in [set('AT'),set('CG')]:counts['palindromic']+=1;continue
  rr=bp[(r.ch,r.pos)];gts={''.join(sorted(x['gt'])) for x in rr if len(x['gt'])==2 and set(x['gt'])<=set('ACGT')}
  if len(gts)!=1:counts['missing_or_ambiguous_call']+=1;continue
  gt=next(iter(gts))
  if not set(gt)<={r.a1,r.a2}:counts['allele_mismatch']+=1;continue
  mapping.append(dict(id=r.id,gt=gt,chrom=r.ch,pos=r.pos,raw_ids=';'.join(x['id'] for x in rr)))
 counts['harmonized']=len(mapping);(A/'hgdp_harmonization.json').write_text(json.dumps(dict(counts=counts,markers=mapping),indent=2))
 (A/'hgdp_extract.txt').write_text('\n'.join(r['id'] for r in mapping)+'\n')
 cmds=[[str(PLINK),'--bfile',str(pref),'--extract',str(A/'hgdp_extract.txt'),'--maf','.05','--geno','.02','--indep-pairwise','200','20','.2','--out',str(A/'hgdp_prune')],[str(PLINK),'--bfile',str(pref),'--extract',str(A/'hgdp_prune.prune.in'),'--make-bed','--out',str(A/'hgdp_reference')]]
 for i,c in enumerate(cmds):
  with (L/f'hgdp_plink_{i}.txt').open('w') as f:subprocess.run(c,stdout=f,stderr=subprocess.STDOUT,check=True)
 print(dict(counts),flush=True)
def run():
 p=A/'hgdp_reference';bim=pd.read_csv(str(p)+'.bim',sep=r'\s+',header=None,names=cols);fam=pd.read_csv(str(p)+'.fam',sep=r'\s+',header=None);meta=pd.read_csv(R/'hgdp650_id_pop_coords.txt',sep='\t').set_index('id').loc[fam[1]]
 n=len(fam);m=len(bim)
 with open(str(p)+'.bed','rb') as f:
  assert f.read(3)==bytes([108,27,1]);buf=np.frombuffer(f.read(),dtype=np.uint8).reshape(m,(n+3)//4)
 codes=((buf[:,:,None]>>np.array([0,2,4,6],dtype=np.uint8))&3).reshape(m,-1)[:,:n];X=np.array([2,np.nan,1,0],dtype=np.float32)[codes].T.copy()
 mapping={r['id']:r for r in json.loads((A/'hgdp_harmonization.json').read_text())['markers']};user=np.array([mapping[r.id]['gt'].count(r.a1) for r in bim.itertuples()],dtype=np.float32)
 labs=meta.population.to_numpy();regions=meta.region.to_numpy();chrom=bim.ch.to_numpy()
 def fit(mask,name,k,markers=None,plot=True):
  xx=X[mask];ls=labs[mask];reg=regions[mask];u=user
  if markers is not None:xx=xx[:,markers];u=u[markers]
  af=np.nanmean(xx,axis=0)/2;ok=(af>.05)&(af<.95);mu=2*af[ok];sd=np.sqrt(2*af[ok]*(1-af[ok]));z=np.nan_to_num((xx[:,ok]-mu)/sd);uz=(u[ok]-mu)/sd
  pca=PCA(n_components=k,svd_solver='randomized',random_state=2026,iterated_power=8);pcs=pca.fit_transform(z);up=pca.transform(uz[None])[0]
  joint=PCA(n_components=k,svd_solver='randomized',random_state=2026,iterated_power=8).fit_transform(np.vstack([z,uz]));jp,jup=joint[:-1],joint[-1]
  def near(pc,uu):
   d=np.linalg.norm(pc-uu,axis=1);idx=np.argsort(d)[:30];return dict(collections.Counter(ls[idx])),[dict(population=ls[i],sample=str(fam[1].to_numpy()[mask][i]),distance=float(d[i])) for i in idx]
  nc,nn=near(pcs,up);jc,jn=near(jp,jup)
  centers=sorted([dict(population=l,distance=float(np.linalg.norm(pcs[ls==l].mean(axis=0)-up)),samples=int(sum(ls==l))) for l in set(ls)],key=lambda r:r['distance'])
  tr,te=train_test_split(np.arange(len(xx)),test_size=.2,stratify=ls,random_state=2026);mu2=np.nanmean(xx[tr][:,ok],axis=0);af2=mu2/2;sd2=np.sqrt(np.maximum(2*af2*(1-af2),.01));zz=np.nan_to_num((xx[:,ok]-mu2)/sd2)
  pp=PCA(n_components=k,svd_solver='randomized',random_state=2026,iterated_power=8);trp=pp.fit_transform(zz[tr]);tep=pp.transform(zz[te]);pred=KNeighborsClassifier(n_neighbors=7).fit(trp,ls[tr]).predict(tep);cv=float(np.mean(pred==ls[te]))
  out=dict(markers=int(ok.sum()),reference_samples=len(xx),populations=len(set(ls)),nearest30_projection=nc,nearest30_joint=jc,nearest_projected_samples=nn,centroids=centers,heldout_population_accuracy=cv,components=k,projection_pc=up.tolist(),joint_pc=jup.tolist())
  pd.DataFrame(pcs,columns=[f'PC{i+1}' for i in range(k)]).assign(population=ls,sample=fam[1].to_numpy()[mask]).to_csv(A/f'{name}_pcs.tsv',sep='\t',index=False)
  if plot:
   fig,ax=plt.subplots(figsize=(10,7),layout='constrained');show=reg if name=='hgdp_global' else ls
   for i,l in enumerate(sorted(set(show))):
    w=show==l;ax.scatter(pcs[w,0],pcs[w,1],s=20,alpha=.65,label=l.replace('_',' '),color=plt.colormaps['tab20'](i))
   ax.scatter(up[0],up[1],s=220,marker='*',color='#d52f54',edgecolor='black',label='Your projected sample',zorder=5)
   ax.set(xlabel='PC1',ylabel='PC2',title=f'{name.replace("_"," ").title()} • {int(ok.sum()):,} SNPs');ax.legend(fontsize=8,loc='center left',bbox_to_anchor=(1,0.5));ax.spines[['top','right']].set_visible(False);fig.savefig(O/f'{name}.png',dpi=150);plt.close(fig)
  print(name,'markers',int(ok.sum()),'projected',nc,'joint',jc,'CV',cv,flush=True);return out
 we=(regions=='Europe')|(regions=='Middle_East_and_Arabia')|(labs=='Mozabite')
 result=dict(global_panel=fit(np.ones(n,dtype=bool),'hgdp_global',10),global_without_reference_outlier=fit(fam[1].to_numpy()!='HGDP00980','hgdp_global_outlier_sensitivity',10,plot=False),west_eurasia=fit(we,'hgdp_west_eurasia',6),west_eurasia_odd=fit(we,'hgdp_west_eurasia_odd',6,chrom%2==1,False),west_eurasia_even=fit(we,'hgdp_west_eurasia_even',6,chrom%2==0,False),limits='Panel similarity, not genealogical percentages. Unequal population sizes and possible relatedness affect nearest-neighbour counts; both centroid and joint-PCA checks retained. Reference-only heldout accuracy is not an individual ancestry probability. HGDP00980 is a reference outlier relative to its Biaka label; missingness 0.000615 does not explain it. Its cause was not established; a separate exclusion sensitivity is retained, not an inference of African ancestry.')
 (A/'hgdp_summary.json').write_text(json.dumps(result,indent=2))
if __name__=='__main__':
 if not (A/'hgdp_reference.bed').exists():prepare()
 run()
