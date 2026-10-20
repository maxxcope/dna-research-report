"""Transparent evidence triage and synthesis, separate from source matching.

This layer does not reclassify variants, compute personal absolute risk, infer
inheritance from two unphased sites, or rate studies using citation counts.
"""
from collections import Counter, defaultdict
import re
from .assay import finding_warning_records

VERSION='1.1'
CLINVAR_REVIEW='https://www.ncbi.nlm.nih.gov/clinvar/docs/review_status/'
ACMG='https://pubmed.ncbi.nlm.nih.gov/25741868/'
CPIC='https://files.cpicpgx.org/data/guideline/publication/serotonin_reuptake_inhibitor_antidepressants/2023/37032427.pdf'
BUCKETS={
 'clinical_review':('Well-reviewed clinical annotation',0),
 'assay_caution':('Assay caveat before interpretation',1),
 'conflicting':('Conflicting clinical evidence',2),
 'limited_review':('Limited clinical review',3),
 'trait':('Curated trait interpretation',4),
 'conditional':('Conditional source pattern',5),
 'medication':('Medication-related source context',6),
 'uncertain':('Uncertain clinical significance',7),
}


def source_warning(finding):
    return ' | '.join(dict.fromkeys(r['warning'] for r in finding_warning_records(finding)))


def make_report(report,genome):
    knowledge=report['knowledge'];findings=knowledge['findings'];by_id={f['id']:f for f in findings}
    warning_records={f['id']:finding_warning_records(f) for f in findings}
    caveats={f['id']:source_warning(f) for f in findings}
    cards=[];groups=defaultdict(list);clinical_ids=set();dedup=set();background=Counter()
    for r in report['clinvar']:
        if warning_records.get(r['rsid']):
            r['source_assay_warnings']=warning_records[r['rsid']]
            r['audit_flags']=sorted(set(r.get('audit_flags',[])+['source_assay_warning']))
        # One allele is not several discoveries because of repeated records.
        signature=(r['chrom'],r['pos'],r['ref'],r['alt'],r['classification'],r['conditions'])
        if signature in dedup:continue
        dedup.add(signature);category=r['category'];clinical_ids.add(r['rsid'])
        flags=r.get('audit_flags',[])
        assay='published_array_quality_flag' in flags or caveats.get(r['rsid'],False)
        if category=='Review with a professional':bucket='assay_caution' if assay else 'clinical_review' if r['stars']>=2 else 'limited_review'
        elif category=='Conflicting evidence':bucket='conflicting'
        elif category=='Uncertain significance':bucket='uncertain'
        else:background[category]+=1;continue
        # Keep multi-gene records together instead of attributing one condition
        # to each gene or inferring an inheritance pattern.
        gene=re.sub(r':\d+','',r['gene']).replace('|',', ') or r['rsid'] or 'Unassigned gene'
        groups[(gene,bucket)].append(r)
    for (gene,bucket),rows in groups.items():
        best=max(r['stars'] for r in rows)
        statements={
          'clinical_review':'The file contains an allele with a pathogenic or likely pathogenic classification supported by an expert panel, guideline, or agreeing submitters. This supports reviewing the specific allele and condition; it does not establish that you have the condition.',
          'assay_caution':'A pathogenic or likely pathogenic source annotation overlaps an assay-quality flag or an explicit source warning about genotyping. Resolve the assay issue before drawing a personal health conclusion.',
          'limited_review':'The source labels an observed allele pathogenic or likely pathogenic, but its review support is limited. Keep it for qualified review rather than treating the label as a diagnosis.',
          'conflicting':'Submitted clinical classifications disagree. This report preserves that disagreement and does not turn it into an increased-risk or disease conclusion.',
          'uncertain':'The clinical meaning is unresolved. This is not evidence that you have the condition, and is not a basis for treatment decisions.',
        }
        variants=[]
        for r in rows:
            f=by_id.get(r['rsid'],{})
            variants.append(dict(id=r['rsid'],variation_id=r['variation_id'],genotype=r['genotype'],alt_copies=r['alt_copies'],
              classification=r['classification'],review=r['review'],review_stars=r['stars'],conditions=r['conditions'],
              source=r['source'],input_markers=r['marker_ids'],flags=r.get('audit_flags',[]),
              source_caution=caveats.get(r['rsid'],'') or '',
              assay_warnings=warning_records.get(r['rsid'],[]),
              context='Condition labels are aggregated; inheritance, penetrance and phase are not inferred.'))
        cards.append(dict(id='clinical:'+gene+':'+bucket,title=gene,bucket=bucket,statement=statements[bucket],
          evidence=f'Highest review among these records: {best}/4 stars. Review stars concern source classification, not personal risk.',
          variants=variants,limits='One allele copy does not by itself establish carrier status or disease; multiple sites do not establish that variants are on different chromosomes.',
          sources=[CLINVAR_REVIEW,ACMG],rank=BUCKETS[bucket][1]))
    for trait in report['traits']:
        if trait['status']!='called':continue
        caution=bool(trait.get('audit_flags'))
        cards.append(dict(id='trait:'+trait['id'],title=trait['title'],bucket='trait',statement=trait['explanation'],
            evidence=trait['evidence'],limits=trait['limits'],sources=trait['sources'],rank=BUCKETS['trait'][1],
            observed=dict(marker=trait['rsid'],genotype=trait['genotype']),assay_caution=caution))
    for rule in knowledge['genoset_audit']:
        partial=rule.get('conditional_evidence')
        if not partial:continue
        cards.append(dict(id='conditional:'+rule['id'],title='Partial evidence: '+rule['summary'],bucket='conditional',
            statement='Some measured markers support this source pattern, but its full rule is not established.',
            evidence=partial['assumption'],limits='This is a conditional source association, not a confirmed trait, risk estimate or medication phenotype.',
            sources=[rule['source']],rank=BUCKETS['conditional'][1],observed_support=partial['observed_support'],
            unobserved_tests=partial['unobserved_tests'],blocking_markers=rule['blocking_markers']))
    medication=defaultdict(list)
    for finding in findings:
        represented=any(v['id']==finding['id'] for card in cards for v in card.get('variants',[]))
        if caveats[finding['id']] and not represented:
            cards.append(dict(id='assay:'+finding['id'],title=(finding.get('gene')+' · ' if finding.get('gene') else '')+finding['id'],bucket='assay_caution',
                statement='A public source flags this input assay or its reported calls as unreliable or error-prone. Keep the source interpretation provisional until the assay issue is resolved.',
                evidence='The warning is retained even when the input probe is represented by a canonical rsID.',
                limits='A probe warning does not establish that a particular measurement is wrong, and a warning page does not establish a genotype or diagnosis.',
                assay_warnings=warning_records[finding['id']],sources=[r['source'] for r in warning_records[finding['id']] if r['source']],rank=BUCKETS['assay_caution'][1]))
    for finding in findings:
        if finding['kind']!='genoset':continue
        text=finding['summary']+' '+finding.get('genotype_text','')
        for gene in set(re.findall(r'\b(?:CYP\d+[A-Z]\d+|DPYD|TPMT|NUDT15|NAT2|G6PD|SLCO1B1|UGT1A1)\b',text)):
            medication[gene].append(finding)
    for gene,rows in medication.items():
        cards.append(dict(id='medication:'+gene,title=gene+' source patterns',bucket='medication',
            statement='The measured markers match the source patterns listed below. These are not a complete medication-metabolism phenotype.',
            evidence='A source haplotype label can describe only part of a gene. Untested alleles, structural variation, copy number and phase can change a complete result.',
            limits='No medication dose or change is inferred. Targeted TPMT phasing, when available, is shown separately with its coverage.',
            source_patterns=[dict(id=r['id'],summary=r['summary'],source=r['source']) for r in rows],
            sources=[CPIC] if gene.startswith('CYP') else [r['source'] for r in rows],rank=BUCKETS['medication'][1]))
    for finding in findings:
        finding['repute']=finding.get('repute') or 'Unknown'
        finding['evidence_scope']='source_multi_marker_rule' if finding['kind']=='genoset' else 'source_association_not_independently_graded'
        finding['source_assay_caution']=bool(caveats[finding['id']])
        finding['source_assay_warning']=caveats[finding['id']]
        finding['assay_warnings']=warning_records[finding['id']]
        finding['clinical_context_available']=finding['id'] in clinical_ids
        matching=[card for card in cards if any(v['id']==finding['id'] for v in card.get('variants',[]))]
        finding['evidence_priority']=min((c['rank'] for c in matching),default=BUCKETS['assay_caution'][1] if caveats[finding['id']] else 20 if finding['kind']=='genoset' else 30)
        finding['interpretation_note']='This is a source pattern, not a complete medication phenotype.' if finding['kind']=='genoset' and re.search(r'metabol|DPYD|TPMT|CYP',finding['summary'],re.I) else 'Source association. Magnitude is an interest rating, not an evidence grade or personal risk.'
        finding['topic']='Medication' if re.search(r'\b(?:metaboliz\w*|metabolis\w*|drug|medication|CYP\w+|DPYD|TPMT)\b',finding['summary']+' '+finding['gene'],re.I) else 'Source associations'
    # A descriptive observation, not an inferred sex, karyotype, or identity.
    y_calls=[r for r in genome.rows if r.chrom=='Y' and set(r.gt)<=set('ACGT') and r.marker not in genome.bad_ids]
    y_evidence=dict(called_markers=len({r.marker for r in y_calls}),meaning='Observed Y-labelled calls are evidence present in this export. No chromosome copy number or sex-development inference is made.')
    cards.sort(key=lambda c:(c['rank'],c['title'],c['id']))
    return dict(version=VERSION,cards=cards,counts=dict(Counter(c['bucket'] for c in cards)),background_annotations=dict(background),
        y_marker_evidence=y_evidence,bucket_labels={k:v[0] for k,v in BUCKETS.items()},
        method='Deterministic evidence triage; source classifications are preserved, not recalculated. No LLM, overall health score, or personal absolute-risk estimate.',
        evidence_limits='SNPedia magnitude and publication counts are not evidence grades. Agreement with SNPedia-based reports is not independent scientific validation.',
        methods_sources=[CLINVAR_REVIEW,ACMG,CPIC])
