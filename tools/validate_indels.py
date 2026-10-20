"""Independent synthetic indel checks against bcftools (via pysam).

Requires the research extra. Writes only an aggregate JSON validation result.
The fixture is generated independently of any genome or saved report.
"""
import argparse
import json
from pathlib import Path
import random
import tempfile
import pysam
import pysam.bcftools
from genetics_app.variants import trim_event


def validate():
    rng=random.Random(1729)
    sequence=''.join(rng.choices('ACGT',k=12000))
    # Repetitive sequence forces bcftools to shift some representations.
    sequence=sequence[:4000]+'AT'*1000+sequence[6000:]
    records=[]
    for i in range(1000):
        anchor=rng.randrange(20,len(sequence)-40);length=rng.randrange(1,15)
        ref=sequence[anchor:anchor+1];alt=ref
        if i%2:ref=sequence[anchor:anchor+length+1]
        else:alt+=('AT'*(length//2+1))[:length] if 4000<anchor<5900 else ''.join(rng.choices('ACGT',k=length))
        # Redundant right padding tests trimming without a shared algorithm.
        if i%3==0:
            tail=sequence[anchor+len(ref):anchor+len(ref)+3];ref+=tail;alt+=tail
        records.append((str(i),anchor+1,ref,alt))
    with tempfile.TemporaryDirectory() as temp:
        root=Path(temp);fa=root/'reference.fa';fa.write_text('>1\n'+sequence+'\n');pysam.faidx(str(fa))
        vcf=root/'input.vcf';vcf.write_text('##fileformat=VCFv4.2\n##contig=<ID=1,length=12000>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n'+''.join(f'1\t{p}\t{i}\t{r}\t{a}\t.\t.\t.\n' for i,p,r,a in sorted(records,key=lambda x:x[1])))
        output=pysam.bcftools.norm('-f',str(fa),'-c','e',str(vcf))
        normalized={}
        for line in output.splitlines():
            if not line.startswith('#'):
                f=line.split('\t');normalized[f[2]]=(int(f[1]),f[3],f[4])
        shifted=0
        for i,pos,ref,alt in records:
            p,deleted,inserted=trim_event(pos,ref,alt)
            assert sequence[p-1:p-1+len(deleted)]==deleted
            application=sequence[:p-1]+inserted+sequence[p-1+len(deleted):]
            q,r,a=normalized[i]
            independent=sequence[:q-1]+a+sequence[q-1+len(r):]
            assert application==independent,(i,pos,ref,alt,q,r,a)
            shifted+=q!=pos
        # An incorrect reference allele must fail the independent check.
        bad=root/'bad.vcf';bad.write_text(vcf.read_text().split('#CHROM')[0]+'#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n1\t1\tbad\t'+('A' if sequence[0]!='A' else 'C')+'\tG\t.\t.\t.\n')
        rejected=False
        try:pysam.bcftools.norm('-f',str(fa),'-c','e',str(bad))
        except pysam.utils.SamtoolsError:rejected=True
        assert rejected
    return dict(fixtures=len(records),equivalent_haplotypes=len(records),shifted_by_bcftools=shifted,wrong_reference_rejected=rejected,
                pysam_version=pysam.__version__,bcftools_version=pysam.__samtools_version__,seed=1729,
                scope='Synthetic allele trimming versus independent bcftools normalization. Does not validate clinical conclusions, source coordinates, or repeat-shifted source allele mapping.')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);args=p.parse_args()
    result=validate();args.out.parent.mkdir(parents=True,exist_ok=True);args.out.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
