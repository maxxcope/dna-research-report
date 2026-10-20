import copy
import hashlib
import json
import sqlite3
import unittest
import test_full_pipeline as fixtures
from genetics_app.assay import warning_passage,probe_warning_records,finding_warning_records
from genetics_app.conclusions import make_report
from genetics_app.curation import apply_definition
from genetics_app.knowledge import interpret
from test_interpretation import report,clinical,genome


class WarningTests(unittest.TestCase):
    def test_detects_probe_wording_and_preserves_clause(self):
        for text in ['Unreliable probe','The provider probe is unreliable (prone to false results).','Error-prone by provider','Often miscalled']:
            self.assertTrue(warning_passage(text))
        pages={'i9001(a;a)':dict(text='',summary='Very prone to miscalling; would otherwise indicate a disease',source='https://example.org/probe')}
        row=probe_warning_records({},pages)['i9001'][0]
        self.assertNotIn('disease',row['warning']);self.assertIn('not asserted to match',row['scope'])
        pages['i9001(a;a)']['summary']='A disease is mentioned, but the probe is unreliable; more disease context'
        row=probe_warning_records({},pages)['i9001'][0]
        self.assertIn('unreliable',row['warning']);self.assertNotIn('disease',row['warning'])

    def test_warning_moves_high_review_annotation_to_assay_caution(self):
        for text in ['Unreliable probe','The provider probe corresponding to this genotype is unreliable (prone to false results).']:
            r=report([clinical(stars=4)]);r['knowledge']['findings'][0]['variant_text']=text
            c=make_report(r,genome('rs9001 1 100 AG\n'))
            self.assertEqual(c['cards'][0]['bucket'],'assay_caution')
            self.assertTrue(c['cards'][0]['variants'][0]['assay_warnings'])
            self.assertTrue(r['clinvar'][0]['source_assay_warnings'])
            self.assertIn('source_assay_warning',r['clinvar'][0]['audit_flags'])

    def test_warning_naming_a_different_probe_does_not_qualify_this_assay(self):
        r=report([clinical()]);f=r['knowledge']['findings'][0]
        f['raw_ids']=['rs9001'];f['variant_text']='The probe i9999 is unreliable.'
        c=make_report(r,genome('rs9001 1 100 AG\n'))
        self.assertEqual(c['cards'][0]['bucket'],'clinical_review')
        self.assertFalse(f['source_assay_caution'])

    def test_general_variant_warning_links_to_its_parent_source_page(self):
        f=dict(id='rs9001',source='https://example.org/rs9001(A;G)',summary='Example',variant_text='Unreliable probe')
        self.assertEqual(finding_warning_records(f)[0]['source'],'https://example.org/rs9001')

    def test_probe_warning_is_visible_without_clinical_annotation(self):
        r=report();r['knowledge']['findings'][0]['probe_assay_warnings']=[dict(source_id='i9001(a;a)',source='https://example.org',warning='Error-prone assay',scope='Probe context only')]
        c=make_report(r,genome('rs9001 1 100 AG\n'))
        self.assertEqual(c['cards'][0]['bucket'],'assay_caution');self.assertEqual(c['cards'][0]['assay_warnings'][0]['source_id'],'i9001(a;a)')
        self.assertEqual(r['knowledge']['findings'][0]['summary'],'Source association')
        self.assertEqual(r['knowledge']['findings'][0]['evidence_priority'],1)


class ProbeMatchingTests(unittest.TestCase):
    setUp=fixtures.KnowledgeTests.setUp
    tearDown=fixtures.KnowledgeTests.tearDown
    genome=fixtures.KnowledgeTests.genome

    def test_alias_warning_survives_without_inheriting_page_genotype(self):
        db=sqlite3.connect(self.db)
        db.execute('insert into genotypes values(?,?)',('i9002(a;a)',json.dumps(dict(summary='Error-prone by provider; would indicate disease',text='',source='https://example.org/probe'))));db.commit();db.close()
        r=interpret(self.genome('i9002 1 200 AG\n'),self.db);f=r['findings'][0]
        self.assertEqual(f['id'],'rs9002');self.assertEqual(f['source_genotype'],'(C;T)')
        self.assertEqual(f['probe_assay_warnings'][0]['warning'],'Error-prone by provider')
        self.assertEqual(f['summary'],'Synthetic minus-strand explanation')
        direct=interpret(self.genome('rs9002 1 200 AG\n'),self.db)['findings'][0]
        self.assertNotIn('probe_assay_warnings',direct)
        no_call=interpret(self.genome('rs9002 1 200 AG\ni9002 1 200 --\n'),self.db)['findings'][0]
        self.assertNotIn('probe_assay_warnings',no_call)


def definition(text,source=('G','A'),forward=('C','T')):
    return dict(id='i9001',gene='SYNTHETIC',name='Synthetic substitution',source_content_sha256=hashlib.sha256(text.encode()).hexdigest(),source_change=list(source),placements={'37':dict(chrom='1',pos=100,ref=forward[0],alt=forward[1])},registry_source='https://example.org/registry',publication='https://example.org/paper')


def marker():return dict(id='i9001',gene='',chrom='',positions={},orientation='',mapping_strand='+',flip37=False,flip38=False)


class CurationTests(unittest.TestCase):
    def test_source_change_and_forward_definition_determine_strand(self):
        text='Public source definition'
        for source,forward,orientation in [(('G','A'),('C','T'),'minus'),(('G','C'),('G','C'),'plus')]:
            m=apply_definition(marker(),text,set(source),definition(text,source,forward))
            self.assertEqual(m['orientations_by_build'],{'37':orientation})
        # Ordered reference/alternate identity resolves palindromic alleles;
        # the person's observed genotype is never an argument to this function.

    def test_changed_source_incompatible_alleles_or_metadata_are_withheld(self):
        text='Public source definition';d=definition(text)
        for m,t,a in [(marker(),'Changed source',{'G','A'}),(marker(),text,{'A'}),({**marker(),'chrom':'2'},text,{'G','A'}),({**marker(),'positions':{'37':101}},text,{'G','A'}),({**marker(),'orientation':'plus'},text,{'G','A'})]:
            self.assertNotIn('curated_snv_definition',apply_definition(m,t,a,d))

    def test_recovered_marker_changes_with_input_and_never_borrows_colocated_no_call(self):
        case=fixtures.KnowledgeTests();case.setUp()
        try:
            db=sqlite3.connect(case.db);m=apply_definition(marker(),'Public source',{'G','A'},definition('Public source'))
            m.update(alias=None,alleles=['A','G'],indel=None,text='Public source',pmids=[],populations=[],retrieved='2026-01-01',source='https://example.org')
            db.execute('insert into markers values(?,?)',('i9001',json.dumps(m)))
            for gt,summary in [('a;a','Alternate genotype'),('a;g','Heterozygote'),('g;g','Misleading condition label')]:
                db.execute('insert into genotypes values(?,?)',('i9001('+gt+')',json.dumps(dict(summary=summary,magnitude=1,repute='Unknown',text='',pmids=[],retrieved='',source='https://example.org'))))
            db.commit();db.close()
            for gt,expected,n in [('CC','(G;G)',0),('CT','(A;G)',1),('TT','(A;A)',2)]:
                f=interpret(genome(f'i9001 1 100 {gt}\n'),case.db)['findings'][0]
                self.assertEqual(f['source_genotype'],expected);self.assertEqual(f['curated_sequence_evidence']['observed_target_alt_copies'],n)
                if not n:self.assertTrue(f['summary'].startswith('Target alternate allele not observed'));self.assertEqual(f['original_source_summary'],'Misleading condition label')
            for rows,reason in [('i9001 1 101 CC\n','coordinate_or_alias_conflict'),('i9001 1 100 AG\n','public_reference_allele_mismatch'),('i9001 1 100 --\nrs9999 1 100 CC\n','no_call')]:
                r=interpret(genome(rows),case.db);self.assertEqual(r['counts']['matched_loci'],0);self.assertEqual(r['unresolved'][0]['reason'],reason)
        finally:case.tearDown()
