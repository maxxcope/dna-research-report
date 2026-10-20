"""Public metadata repair must not use personal alleles to choose a strand."""
import copy
import json
import sqlite3
import unittest
import test_full_pipeline as fixtures
from genetics_app.source_enrichment import unique_snv,enrich_marker,response_records,validate_cache,API,SELECTION_VERSION
from genetics_app.knowledge import import_knowledge,interpret
from genetics_app.references import digest


def record(ch='1',pos=100,alleles='A/G'):
    return dict(name='rs9001',var_class='SNP',synonyms=['rs9991'],mappings=[dict(assembly_name='GRCh37',seq_region_name=ch,start=pos,end=pos,strand=1,allele_string=alleles)])


def marker(source=()):
    return dict(id='rs9001',chrom='1',positions={'37':90,'38':80},orientation='',alleles=list(source),indel=None)


class PublicPlacementTests(unittest.TestCase):
    def test_primary_mapping_distinguishes_alt_scaffolds_and_multiple_primary_sites(self):
        r=record();r['mappings'].append({**r['mappings'][0],'seq_region_name':'ALT_PATCH'})
        self.assertEqual(unique_snv(r)['pos'],100)
        r['mappings'].append({**r['mappings'][0],'seq_region_name':'2'})
        self.assertIsNone(unique_snv(r))

    def test_complex_and_reverse_records_withheld(self):
        for r in [record(alleles='A/AT'),{**record(),'var_class':'indel'}]:self.assertIsNone(unique_snv(r))
        r=record();r['mappings'][0]['strand']=-1;self.assertIsNone(unique_snv(r))

    def test_multiallelic_SNV_placement_does_not_imply_strand(self):
        p=unique_snv(record(alleles='A/C/G'));self.assertEqual(p['pos'],100)
        m=enrich_marker(marker(),p,{'A','G'});self.assertEqual(m['orientations_by_build'],{'37':'plus'})
        p=unique_snv(record(alleles='A/C/G/T'))
        self.assertNotIn('orientations_by_build',enrich_marker(marker(),p,{'A','G'}))
        self.assertNotIn('orientations_by_build',enrich_marker(marker(),p,{'A','T'}))

    def test_only_documented_response_aliases_are_accepted(self):
        r=record();self.assertEqual(response_records({'rs9001':r},['rs9991']),{'rs9991':r})
        with self.assertRaises(ValueError):response_records({'rs9001':r},['rs9992'])

    def test_orientation_derived_from_full_source_alleles(self):
        p=unique_snv(record())
        for alleles,expected in [({'A','G'},'plus'),({'T','C'},'minus')]:
            m=enrich_marker(marker(),p,alleles)
            self.assertEqual(m['orientations_by_build'],{'37':expected});self.assertEqual(m['positions'],{'37':100,'38':80})
        for alleles in [set(),{'A'},{'A','C','G','T'}]:self.assertNotIn('orientations_by_build',enrich_marker(marker(),p,alleles))
        self.assertNotIn('orientations_by_build',enrich_marker(marker(),unique_snv(record(alleles='A/T')),{'A','T'}))

    def test_conflicting_chromosome_and_indel_definitions_never_overwritten(self):
        m=marker(['-','AG']);self.assertEqual(enrich_marker(copy.deepcopy(m),unique_snv(record()),set()),m)
        m=marker();self.assertEqual(enrich_marker(copy.deepcopy(m),unique_snv(record(ch='2')),{'A','G'}),m)


class EnrichedMatchingTests(unittest.TestCase):
    setUp=fixtures.KnowledgeTests.setUp
    tearDown=fixtures.KnowledgeTests.tearDown
    genome=fixtures.KnowledgeTests.genome

    def prepare(self,r,source_alleles=('A','G'),orientation=''):
        source=self.root/'source.db';db=sqlite3.connect(source)
        db.execute("update snps set content=? where rsid='rs9001'",('{{Rsnum|rsid=9001|Chromosome=1|position=90|Assembly=GRCh38|Orientation='+orientation+'}}',))
        db.execute("delete from genotypes where id like 'rs9001(%'")
        for a in source_alleles:db.execute('insert into genotypes values(?,?,?)',(f'rs9001({a.lower()};{a.lower()})','{{Genotype|summary=Synthetic allele definition}}','2026-01-01'))
        db.commit();db.close()
        cache=self.root/'placements.json';cache.write_text(json.dumps(dict(schema=1,kind='Ensembl public placements',source=API,selection_version=SELECTION_VERSION,complete=True,build='37',database_sha256=digest(source),criteria_sha256=digest(self.root/'criteria.json'),chain_sha256='',requested=['rs9001'],missing=[],records={'rs9001':r},retrieved_at='2026-01-01')))
        import_knowledge(source,self.root/'criteria.json',None,self.db,placement_reference=cache)

    def test_incomplete_mismatched_or_unrelated_reference_cache_rejected(self):
        self.prepare(record())
        state=json.loads((self.root/'placements.json').read_text())
        for change in [dict(complete=False),dict(criteria_sha256='wrong'),dict(records={}),dict(missing=['rs9001']),dict(records={'rs9001':{**record(),'name':'rs9002','synonyms':[]}})]:
            with self.assertRaises(ValueError):validate_cache({**state,**change},self.root/'source.db',self.root/'criteria.json',None)

    def test_orphan_genotype_page_uses_public_placement_without_invented_strand(self):
        self.prepare(record())
        source=self.root/'source.db';db=sqlite3.connect(source);db.execute("delete from snps where rsid='rs9001'");db.commit();db.close()
        cache=self.root/'placements.json';state=json.loads(cache.read_text());state['database_sha256']=digest(source);cache.write_text(json.dumps(state))
        import_knowledge(source,self.root/'criteria.json',None,self.db,placement_reference=cache)
        result=interpret(self.genome('rs9001 1 100 AG\n'),self.db)
        self.assertEqual(result['findings'][0]['source_genotype'],'(A;G)')

    def test_repaired_placement_and_strand_all_input_genotypes(self):
        self.prepare(record(alleles='C/T'),('A','G'))
        for gt,expected in [('CC','(G;G)'),('CT','(A;G)'),('TT','(A;A)')]:
            f=interpret(self.genome(f'rs9001 1 100 {gt}\n'),self.db)['findings'][0]
            self.assertEqual(f['source_genotype'],expected)
        self.assertEqual(interpret(self.genome('rs9001 1 101 CT\n'),self.db)['counts']['matched_loci'],0)
        self.assertEqual(interpret(self.genome('rs9001 1 100 AG\n'),self.db)['unresolved'][0]['reason'],'public_reference_allele_mismatch')

    def test_palindromic_heterozygote_is_invariant_but_homozygote_stays_unknown(self):
        self.prepare(record(alleles='A/T'),('A','T'))
        result=interpret(self.genome('rs9001 1 100 AT\n'),self.db)
        self.assertEqual(result['findings'][0]['source_genotype'],'(A;T)')
        self.assertIn('strand_invariant_heterozygote',result['findings'][0]['flags'])
        for gt in ('AA','TT'):self.assertEqual(interpret(self.genome(f'rs9001 1 100 {gt}\n'),self.db)['counts']['matched_loci'],0)

    def test_explicit_public_merged_id_can_supply_rule_call(self):
        self.prepare(record())
        f=interpret(self.genome('rs9991 1 100 AG\n'),self.db)['findings'][0]
        self.assertEqual(f['id'],'rs9001');self.assertEqual(f['public_aliases'],['rs9991'])
        self.assertEqual(interpret(self.genome('rs9992 1 100 AG\n'),self.db)['counts']['matched_loci'],0)
        self.assertEqual(interpret(self.genome('rs9001 1 100 AA\nrs9991 1 100 GG\n'),self.db)['counts']['matched_loci'],0)

    def test_alias_with_incompatible_coordinate_cannot_supply_or_override_a_call(self):
        self.prepare(record())
        r=interpret(self.genome('rs9991 1 101 AG\n'),self.db)
        self.assertEqual(r['counts']['matched_loci'],0)
        self.assertEqual(r['unresolved'][0]['reason'],'documented_alias_coordinate_conflict')
        r=interpret(self.genome('rs9001 1 100 AA\nrs9991 1 101 GG\n'),self.db)
        self.assertEqual(r['findings'][0]['source_genotype'],'(A;A)')
        self.assertEqual(r['findings'][0]['excluded_public_aliases'],['rs9991'])
