"""Public, allele-specific normalization for simple insertion/deletion assays.

Only explicit rsID/alias identities and published coordinates are accepted.
No nearest-position search, generic +/-1 tolerance, or I/D sequence guessing.
"""
from collections import defaultdict
import re
from .references import connect_read


def encoding_warning(marker):
    # A published exception to the provider's allele coding defeats generic
    # I=long/D=short conversion. Preserve the exception, never guess a reversal.
    return bool(re.search(r'(?:incorrectly|swapp?ed|reversed).{0,100}(?:assign|genotype|allele|deletion|insertion)',marker.get('text',''),re.I))


def reverse_complement(sequence):
    return sequence.translate(str.maketrans('ACGT','TGCA'))[::-1]


def trim_event(pos, ref, alt):
    """Return the unpadded 1-based event; not repeat left-normalization."""
    if not ref or not alt or not set(ref+alt)<=set('ACGT'):return None
    while ref and alt and ref[0]==alt[0]:
        pos+=1;ref=ref[1:];alt=alt[1:]
    while ref and alt and ref[-1]==alt[-1]:ref=ref[:-1];alt=alt[:-1]
    if bool(ref)==bool(alt):return None  # Exclude substitutions/delins/no change.
    return pos,ref,alt


def public_indels(path,wanted):
    """Unique public alleles, including clinical context, independently of input."""
    if not path:return {},None
    db,meta=connect_read(path);records=defaultdict(dict);rejected=set()
    try:
        if meta.get('kind')!='ClinVar':raise ValueError('Expected public ClinVar catalogue.')
        for row in db.execute('select * from variants'):
            r=dict(row)
            for token in re.split('[;,]',r['rs']):
                key='rs'+token
                if key not in wanted:continue
                event=trim_event(r['pos'],r['ref'],r['alt'])
                if event is None:rejected.add(key);continue
                signature=(r['chrom'],r['pos'],r['ref'],r['alt'])
                records[key][signature]={**r,'event_pos':event[0],'deleted':event[1],'inserted':event[2]}
    finally:db.close()
    return {key:next(iter(rows.values())) for key,rows in records.items() if len(rows)==1 and key not in rejected},meta['build']


def source_alleles(fields):
    alleles=set()
    for key in ('geno1','geno2','geno3'):
        match=re.fullmatch(r'\(([ACGT-]+);([ACGT-]+)\)',fields.get(key,'').upper())
        if match:alleles.update(match.groups())
    return sorted(alleles)


def prepare_indel(marker,record,build):
    if not record or marker['orientation'] not in {'plus','minus'}:return None
    alleles=marker.get('alleles',[])
    if len(alleles)!=2 or '-' not in alleles:return None
    long=next(a for a in alleles if a!='-')
    if not set(long)<=set('ACGT'):return None
    forward=reverse_complement(long) if marker['orientation']=='minus' else long
    if {record['deleted'],record['inserted']}!={'',forward}:return None
    if marker['chrom']!=record['chrom'] or marker['mapping_strand']=='-' or marker['flip'+build]:return None
    # These are explicit source/VCF/event conventions for this exact allele,
    # not arbitrary neighbouring bases.
    positions={record['pos'],record['event_pos']}
    if marker['positions'].get(build):positions.add(marker['positions'][build])
    return dict(build=build,positions=sorted(positions),source_insertion=long,forward_insertion=forward,
                method='unique_ClinVar_allele_and_SNPedia_sequence_agreement',clinical_record=record)


def resolve_indel(genome,marker,ids):
    definition=marker.get('indel')
    if encoding_warning(marker):return None
    if not definition or definition['build']!=genome.build or genome.qc.get('provider')!='23andMe':return None
    rows=[r for key in ids for r in genome.by_id.get(key,[])]
    if not rows or any(key in genome.bad_ids for key in ids):return None
    if any(r.chrom!=marker['chrom'] or r.pos not in definition['positions'] for r in rows):return None
    calls={r.gt for r in rows if r.gt!='--'}
    if len(calls)!=1:return None
    gt=next(iter(calls))
    if not set(gt)<=set('ID') or (marker['chrom'].isdigit() and len(gt)!=2):return None
    source=tuple(sorted(definition['source_insertion'] if a=='I' else '-' for a in gt))
    sequence=tuple(definition['forward_insertion'] if a=='I' else '' for a in gt)
    return dict(raw_genotype=gt,source_alleles=source if len(source)==2 else source*2,
                sequence_call=sequence,raw_positions=sorted({r.pos for r in rows}),definition=definition,
                coding_status='polarity_invariant_heterozygote' if set(gt)=={'I','D'} else 'provider_coding_unverified')
