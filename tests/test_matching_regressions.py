"""Counterexamples for reference rescue, assay identity, and missing evidence."""
import json
import sqlite3
import unittest
import test_full_pipeline as fixtures
from genetics_app.knowledge import import_knowledge, interpret, public_snv_positions
from genetics_app.rules import evaluate_all, unresolved_markers


class MatchingRegressions(unittest.TestCase):
    setUp=fixtures.KnowledgeTests.setUp
    tearDown=fixtures.KnowledgeTests.tearDown
    genome=fixtures.KnowledgeTests.genome

    def test_colocated_indel_does_not_invalidate_snv(self):
        result=interpret(self.genome('rs9001 1 100 AG\ni9999 1 100 II\n'),self.db)
        row=result['findings'][0]
        self.assertEqual(row['source_genotype'],'(A;G)')
        self.assertEqual(row['raw_ids'],['rs9001'])
        self.assertEqual(row['public_aliases'],[])

    def test_colocated_disagreeing_snv_still_blocks(self):
        result=interpret(self.genome('rs9001 1 100 AG\ni9999 1 100 II\nrs9998 1 100 AA\n'),self.db)
        self.assertEqual(result['counts']['matched_loci'],0)

    def test_unrelated_colocated_probe_cannot_rescue_no_call(self):
        result=interpret(self.genome('rs9001 1 100 --\nrs9999 1 100 AG\n'),self.db)
        self.assertEqual(result['counts']['matched_loci'],0)
        self.assertEqual(result['unresolved'][0]['reason'],'no_call')

    def test_official_mitochondrial_positions_are_identical_across_builds(self):
        import gzip
        source=self.root/'source.db';db=sqlite3.connect(source)
        db.execute("update snps set content='{{Rsnum|Chromosome=MT|position=100|Assembly=GRCh38|Orientation=plus}}' where rsid='rs9001'")
        db.commit();db.close()
        chain=self.root/'chain.gz'
        with gzip.open(chain,'wt') as f:f.write('chain 1 chrMT 200 + 0 100 chrMT 200 + 10 110 1\n100\n')
        import_knowledge(source,self.root/'criteria.json',chain,self.db)
        db=sqlite3.connect(self.db);m=json.loads(db.execute("select data from markers where id='rs9001'").fetchone()[0]);db.close()
        self.assertEqual(m['positions'],{'37':100,'38':100})

    def test_public_indel_rs_with_another_allele_is_rejected(self):
        from genetics_app.variants import public_indels
        path=self.root/'indel.db';db=sqlite3.connect(path)
        db.executescript("create table metadata(key text,value text);insert into metadata values('schema','1'),('kind','ClinVar'),('build','37');create table variants(chrom text,pos integer,ref text,alt text,rs text);")
        db.executemany('insert into variants values(?,?,?,?,?)',[
            ('1',100,'GA','G','1'),('1',100,'G','A','1'),
            ('1',200,'GA','G','2'),('1',201,'AT','A','2'),('1',300,'GA','G','3')]);db.commit();db.close()
        self.assertEqual(set(public_indels(path,{'rs1','rs2','rs3'})[0]),{'rs3'})

    def test_indel_needs_literal_source_page_and_provider(self):
        from genetics_app.parser import parse_genome
        db=sqlite3.connect(self.db)
        row=json.loads(db.execute("select data from genotypes where id='rs9001(a;g)'").fetchone()[0])
        row['id']='rs9001(i;i)';db.execute('insert into genotypes values(?,?)',(row['id'],json.dumps(row)));db.commit();db.close()
        yes=interpret(parse_genome(b'# 23andMe build 37\nrs9001 1 100 II\n'),self.db)
        self.assertEqual(yes['findings'][0]['source_genotype'],'(I;I)')
        for data in [b'# build 37\nrs9001 1 100 II\n',b'# 23andMe build 37\nrs9001 1 100 DD\n',b'# 23andMe build 37\nrs9001 1 100 I\n']:
            self.assertEqual(interpret(parse_genome(data),self.db)['counts']['matched_loci'],0)

    def test_unverified_homozygous_code_cannot_generate_sequence_finding(self):
        from test_interpretation import marker,genome
        db=sqlite3.connect(self.db)
        original=json.loads(db.execute("select data from markers where id='rs9001'").fetchone()[0])
        m={**original,**marker()};db.execute("update markers set data=? where id='rs9001'",(json.dumps(m),));db.commit();db.close()
        for gt in ('II','DD'):
            result=interpret(genome(f'rs9001 7 101 {gt}\n'),self.db)
            self.assertEqual(result['counts']['matched_loci'],0)
            self.assertEqual(result['unresolved'][0]['reason'],'provider_indel_coding_unverified')
        result=interpret(genome('rs9001 7 101 DI\n'),self.db)
        self.assertEqual(result['findings'][0]['source_genotype'],'(-;A)')

    def test_public_coordinate_rescue_not_nearest_position(self):
        path=self.root/'clinvar.db';db=sqlite3.connect(path)
        db.executescript("create table metadata(key text,value text);insert into metadata values('schema','1'),('kind','ClinVar'),('build','37');create table variants(chrom text,pos integer,ref text,alt text,rs text);")
        db.executemany('insert into variants values(?,?,?,?,?)',[
            ('1',101,'A','G','9001'),('1',201,'A','C','9002'),('1',202,'A','C','9002'),
            ('1',301,'A','AT','9003')]);db.commit();db.close()
        locations,build=public_snv_positions(path,{'rs9001','rs9002','rs9003'})
        self.assertEqual(locations,{'rs9001':('1',101)})
        import_knowledge(self.root/'source.db',self.root/'criteria.json',None,self.db,path)
        yes=interpret(self.genome('rs9001 1 101 AG\n'),self.db)
        self.assertIn('source_coordinate_corrected_from_ClinVar',yes['findings'][0]['flags'])
        self.assertEqual(interpret(self.genome('rs9001 1 102 AG\n'),self.db)['counts']['matched_loci'],0)

    def test_unknown_negative_is_not_healthy_and_trace_excludes_decided_branches(self):
        rules={'gs1':'and(rs1(A),not(rs2(T)),or(rs3(G),rs4(C)))'}
        calls={'rs1':('A','A'),'rs3':('G','G')}
        outcomes,_=evaluate_all(rules,calls)
        self.assertIsNone(outcomes['gs1'])
        self.assertEqual(unresolved_markers('gs1',rules,calls,outcomes),['rs2'])
        calls['rs2']=('C','C');self.assertTrue(evaluate_all(rules,calls)[0]['gs1'])
        calls['rs2']=('T','T');self.assertFalse(evaluate_all(rules,calls)[0]['gs1'])


if __name__=='__main__':unittest.main()
