import json
from pathlib import Path
import tempfile
import unittest
from genetics_app.rules import evaluate_all, conditional_evidence
from genetics_app.variants import trim_event, prepare_indel, resolve_indel
from genetics_app.parser import parse_genome
from genetics_app.conclusions import make_report
from genetics_app.engine import classification_category, sequence_clinvar_results


def marker(orientation='plus',seq='A'):
    m=dict(id='rs9001',chrom='7',positions={'37':102},orientation=orientation,mapping_strand='+',flip37=False,flip38=False,alleles=['-',seq])
    row=dict(id='9001',chrom='7',pos=100,ref='GA',alt='G',event_pos=101,deleted='A',inserted='',rs='9001',gene='TEST:1',classification='Pathogenic',review='reviewed_by_expert_panel',conditions='Example',origin='1',conflicts='')
    m['indel']=prepare_indel(m,row,'37');return m


def genome(rows,provider=True):
    return parse_genome(('# '+('23andMe ' if provider else '')+'build 37\n'+rows).encode())


class VariantTests(unittest.TestCase):
    def test_vcf_padding_and_complex_rejection(self):
        self.assertEqual(trim_event(100,'GA','G'),(101,'A',''))
        self.assertEqual(trim_event(100,'G','GAT'),(101,'','AT'))
        self.assertIsNone(trim_event(100,'ATCTGA','TAGATTC'))
        self.assertIsNone(trim_event(100,'A','<DEL>'))
        self.assertIsNone(trim_event(100,'A','C'))

    def test_only_explicit_public_positions(self):
        m=marker()
        for p in (100,101,102):
            r=resolve_indel(genome(f'rs9001 7 {p} DI\n'),m,{'rs9001'})
            self.assertEqual(r['source_alleles'],('-','A'))
        for p in (99,103):self.assertIsNone(resolve_indel(genome(f'rs9001 7 {p} DI\n'),m,{'rs9001'}))

    def test_reject_uncorroborated_sequence_provider_conflicts_and_ploidy(self):
        self.assertIsNone(marker(seq='T')['indel'])
        self.assertIsNone(resolve_indel(genome('rs9001 7 101 II\n',False),marker(),{'rs9001'}))
        self.assertIsNone(resolve_indel(genome('rs9001 7 101 I\n'),marker(),{'rs9001'}))
        self.assertIsNone(resolve_indel(genome('rs9001 7 101 II\ni9001 7 101 DD\n'),marker(),{'rs9001','i9001'}))

    def test_reverse_complement_uses_whole_sequence(self):
        from genetics_app.variants import reverse_complement
        self.assertEqual(reverse_complement('AGTC'),'GACT')
        m=marker('minus','T')
        r=resolve_indel(genome('rs9001 7 101 DI\n'),m,{'rs9001'})
        self.assertEqual(r['source_alleles'],('-','T'));self.assertEqual(r['sequence_call'],('','A'))

    def test_published_coding_exception_withholds_conversion(self):
        m=marker();m['text']='The provider incorrectly assigns the deletion genotype as the normal genotype.'
        self.assertIsNone(resolve_indel(genome('rs9001 7 101 DD\n'),m,{'rs9001'}))

    def test_clinvar_alt_copy_semantics_deletion(self):
        for gt,expected in [('II',0),('DI',1),('DD',0)]:
            r=resolve_indel(genome(f'rs9001 7 101 {gt}\n'),marker(),{'rs9001'})
            rows=sequence_clinvar_results({'findings':[dict(id='rs9001',raw_ids=['rs9001'],flags=[],sequence_resolution=r)]})
            if expected:self.assertEqual(rows[0]['alt_copies'],expected)
            else:self.assertEqual(rows,[])

    def test_insertion_alt_is_I_not_D(self):
        m=marker();record=m['indel']['clinical_record'];record.update(ref='G',alt='GA',deleted='',inserted='A');m['indel']=prepare_indel(m,record,'37')
        for gt,expected in [('II',0),('DI',1),('DD',0)]:
            r=resolve_indel(genome(f'rs9001 7 101 {gt}\n'),m,{'rs9001'})
            rows=sequence_clinvar_results({'findings':[dict(id='rs9001',raw_ids=['rs9001'],flags=[],sequence_resolution=r)]})
            self.assertEqual(rows[0]['alt_copies'] if rows else 0,expected)

    def test_both_strands_and_all_diploid_code_combinations(self):
        for orientation,source in [('plus','A'),('minus','T')]:
            for gt in ['II','ID','DI','DD']:
                r=resolve_indel(genome(f'rs9001 7 101 {gt}\n'),marker(orientation,source),{'rs9001'})
                self.assertEqual(r['source_alleles'].count('-'),gt.count('D'))
                self.assertEqual(r['source_alleles'].count(source),gt.count('I'))
                self.assertEqual(r['coding_status'],'polarity_invariant_heterozygote' if len(set(gt))==2 else 'provider_coding_unverified')


class ConditionalTests(unittest.TestCase):
    def test_missing_exclusion_is_explained_not_matched(self):
        rules={'gs1':'and(rs1(A),not(rs2(T)))'};calls={'rs1':('A','A')}
        strict,_=evaluate_all(rules,calls);self.assertIsNone(strict['gs1'])
        partial=conditional_evidence(rules,calls,strict)['gs1']
        self.assertEqual(partial['unobserved_tests'],['rs2(T)'])
        self.assertEqual(partial['observed_support'][0]['id'],'rs1')
        for gt,expected in [(('T','T'),False),(('C','C'),True)]:
            complete={**calls,'rs2':gt};out,_=evaluate_all(rules,complete)
            self.assertEqual(out['gs1'],expected);self.assertEqual(conditional_evidence(rules,complete,out),{})

    def test_no_false_normal_from_absent_data(self):
        for rule in ['not(rs2(T))','and(rs1(A),not(rs2(T)))','or(not(rs2(T)),not(rs3(A)))']:
            rules={'gs1':rule};out,_=evaluate_all(rules,{})
            self.assertEqual(conditional_evidence(rules,{},out),{})

    def test_nested_negative_and_cycle(self):
        rules={'gs1':'and(rs1(A),not(gs2))','gs2':'or(rs2(T),rs3(G))'};calls={'rs1':('A','A'),'rs2':('C','C')}
        out,_=evaluate_all(rules,calls);partial=conditional_evidence(rules,calls,out)
        self.assertEqual(partial['gs1']['unobserved_tests'],['rs3(G)'])
        rules={'gs1':'and(rs1(A),not(gs2))','gs2':'gs1'};out,_=evaluate_all(rules,calls)
        self.assertEqual(conditional_evidence(rules,calls,out),{})


def clinical(category='Review with a professional',stars=3,flags=None):
    return dict(chrom='1',pos=100,ref='A',alt='G',classification='Pathogenic',conditions='Example condition',category=category,stars=stars,
                gene='EXAMPLE:1',rsid='rs9001',genotype='AG',alt_copies=1,variation_id='9001',review='expert panel',source='https://example.org',marker_ids=['rs9001'],audit_flags=flags or [])


def report(rows=None,source_summary='Source association',mag=10):
    f=dict(id='rs9001',kind='variant',summary=source_summary,genotype_text='',variant_text='',repute='Bad',gene='EXAMPLE',flags=[],magnitude=mag)
    return dict(knowledge=dict(findings=[f],genoset_audit=[]),clinvar=rows or [],traits=[])


class ConclusionTests(unittest.TestCase):
    def synthesize(self,r):return make_report(r,genome('rs9001 1 100 AG\n'))
    def test_source_magnitude_never_becomes_clinical_evidence(self):
        self.assertEqual(self.synthesize(report())['cards'],[])
    def test_review_uses_actual_review_not_popularity(self):
        self.assertEqual(self.synthesize(report([clinical(stars=3)]))['cards'][0]['bucket'],'clinical_review')
        self.assertEqual(self.synthesize(report([clinical(stars=0)]))['cards'][0]['bucket'],'limited_review')
    def test_uncertainty_and_conflicts_never_promoted(self):
        for category,bucket in [('Conflicting evidence','conflicting'),('Uncertain significance','uncertain')]:
            r=clinical(category=category,stars=4)
            self.assertEqual(self.synthesize(report([r]))['cards'][0]['bucket'],bucket)
        self.assertEqual(classification_category('Pathogenic;Uncertain_significance'),'Conflicting evidence')
        self.assertEqual(classification_category('Pathogenic; Uncertain significance'),'Conflicting evidence')
        self.assertEqual(classification_category('Pathogenic/Benign'),'Conflicting evidence')
    def test_assay_warning_survives_high_source_review(self):
        self.assertEqual(self.synthesize(report([clinical(stars=4,flags=['published_array_quality_flag'])]))['cards'][0]['bucket'],'assay_caution')
        self.assertEqual(self.synthesize(report([clinical()],source_summary='Likely miscall by platform'))['cards'][0]['bucket'],'assay_caution')
    def test_duplicate_records_not_multiple_findings(self):
        r=clinical();cards=self.synthesize(report([r,r]))['cards'];self.assertEqual(len(cards[0]['variants']),1)

    def test_warning_passage_preserved_from_general_source(self):
        r=report([clinical()],source_summary='Unrelated summary')
        r['knowledge']['findings'][0]['variant_text']='Variant description. This variant is often a miscall. Other details.'
        card=self.synthesize(r)['cards'][0]
        self.assertIn('often a miscall',card['variants'][0]['source_caution'])
        self.assertNotIn('Unrelated summary',card['variants'][0]['source_caution'])
    def test_two_unphased_sites_do_not_imply_compound_heterozygosity(self):
        a=clinical();b={**clinical(),'pos':200,'rsid':'rs9002'};card=self.synthesize(report([a,b]))['cards'][0]
        self.assertEqual(len(card['variants']),2);self.assertIn('do not establish',card['limits'])
    def test_benign_does_not_create_global_all_clear(self):
        r=self.synthesize(report([clinical(category='Benign / likely benign')]))
        self.assertEqual(r['cards'],[]);self.assertEqual(r['background_annotations']['Benign / likely benign'],1)


if __name__=='__main__':unittest.main()
