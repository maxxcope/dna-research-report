import gzip
import hashlib
import json
from pathlib import Path
import sqlite3
import tempfile
import unittest
from unittest.mock import patch
from genetics_app.knowledge import Chain, import_knowledge, interpret, parameters, population_rows
from genetics_app.parser import parse_genome
from genetics_app.scoring import direct_scores
from genetics_app.extras import quality

class KnowledgeTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.root=Path(self.tmp.name)
        source=self.root/'source.db';db=sqlite3.connect(source)
        db.executescript('create table snps(rsid text,content text,scraped_at text);create table genotypes(id text,content text,scraped_at text);create table genosets(id text,content text,scraped_at text);')
        for marker,orientation,pos in [('rs9001','plus',100),('rs9002','minus',200),('rs9003','plus',300)]:
            text='{{Rsnum|rsid='+marker[2:]+'|Chromosome=1|position='+str(pos)+'|Assembly=GRCh37|Orientation='+orientation+'|Gene=SYNTHETIC}}'
            db.execute('insert into snps values(?,?,?)',(marker,text,'2026-01-01'))
        db.execute('insert into snps values(?,?,?)',('i9002','{{23andMe SNP|rsid=9002}}','2026-01-01'))
        db.execute('insert into genotypes values(?,?,?)',('rs9001(a;g)','{{Genotype|summary=Synthetic heterozygous explanation|magnitude=2|repute=Good}}','2026-01-01'))
        db.execute('insert into genotypes values(?,?,?)',('rs9002(c;t)','{{Genotype|summary=Synthetic minus-strand explanation|magnitude=3|repute=Bad}}','2026-01-01'))
        db.execute('insert into genosets values(?,?,?)',('gs9001','{{Genoset|summary=Synthetic combination|Magnitude=4|Repute=Bad}}','2026-01-01'));db.commit();db.close()
        rule='and(rs9001(A;G),rs9002(C;T))';criteria=self.root/'criteria.json';criteria.write_text(json.dumps(dict(complete=True,pages=[dict(title='Gs9001/criteria',revisions=[dict(revid=1,timestamp='2026-01-01',sha1=hashlib.sha1(rule.encode()).hexdigest(),slots={'main':{'content':rule}})])])))
        self.db=self.root/'knowledge.db';import_knowledge(source,criteria,None,self.db)
    def tearDown(self):self.tmp.cleanup()
    def genome(self,rows):return parse_genome(('# build 37\n'+rows).encode())
    def test_same_library_changes_with_input(self):
        yes=interpret(self.genome('rs9001 1 100 AG\nrs9002 1 200 AG\n'),self.db)
        no=interpret(self.genome('rs9001 1 100 AA\nrs9002 1 200 AG\n'),self.db)
        self.assertEqual(yes['counts']['genosets'],1);self.assertEqual(no['counts']['genosets'],0)
        self.assertEqual(yes['counts']['variant_summaries'],2);self.assertEqual(no['counts']['variant_summaries'],1)
    def test_missing_call_not_negative_rule(self):
        result=interpret(self.genome('rs9001 1 100 AG\nrs9002 1 200 --\n'),self.db)
        self.assertEqual(result['genoset_audit'][0]['status'],'unresolved')
    def test_coordinate_mismatch_withheld(self):
        result=interpret(self.genome('rs9001 1 101 AG\nrs9002 1 200 AG\n'),self.db)
        self.assertEqual(result['counts']['genosets'],0)
        self.assertIn('coordinate_or_alias_conflict',[r['reason'] for r in result['unresolved']])
    def test_public_alias_and_quality_propagation(self):
        result=interpret(self.genome('rs9001 1 100 AG\ni9002 1 200 AG\n'),self.db,{('1',200)})
        self.assertEqual(result['counts']['genosets'],1)
        gs=next(r for r in result['findings'] if r['kind']=='genoset');self.assertIn('published_array_quality_flag',gs['flags'])
    def test_merged_probe_conflict_withheld(self):
        result=interpret(self.genome('rs9001 1 100 AG\nrs9002 1 200 AG\ni9002 1 200 AA\n'),self.db)
        self.assertEqual(result['counts']['genosets'],0)
    def test_inline_nested_template_and_population(self):
        text='{{Genotype|summary=A [[Thing|label]] and {{small|nested}}|repute=Good}}'
        self.assertEqual(parameters(text,{'genotype'})['summary'],'A [[Thing|label]] and {{small|nested}}')
        p=population_rows('{{ population diversity|geno1=(A;A)|geno2=(A;G)|geno3=(G;G)|CEU|10|20|70|HapMapRevision=28}}')
        self.assertEqual(p[0]['percentages'],[10,20,70])
    def test_source_code_is_only_text(self):
        db=sqlite3.connect(self.db);r=json.loads(db.execute("select data from genotypes where id='rs9001(a;g)'").fetchone()[0]);r['summary']='<script>alert(1)</script>';db.execute("update genotypes set data=? where id='rs9001(a;g)'",(json.dumps(r),));db.commit();db.close()
        from genetics_app.report import report_lines
        from genetics_app.engine import analyze
        data=b'# build 37\nrs9001 1 100 AG\n';report=analyze(data,traits_only=True);report['knowledge']=interpret(parse_genome(data),self.db)
        text='\n'.join(report_lines(report));self.assertIn('# DNA Research Report',text);self.assertIn('<script>alert(1)</script>',text)

class ChainTests(unittest.TestCase):
    def test_forward_reverse_gaps_and_ambiguity(self):
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/'test.gz'
            with gzip.open(path,'wt') as f:f.write('chain 10 chr1 1000 + 100 110 chr1 2000 + 200 210 1\n10\n\nchain 10 chr2 1000 + 100 110 chr2 2000 - 200 210 2\n10\n')
            c=Chain(path);self.assertEqual(c.map('1',101),(201,'+'));self.assertEqual(c.map('1',201,True),(101,'+'));self.assertIsNone(c.map('1',111));self.assertEqual(c.map('2',101),(1800,'-'));self.assertEqual(c.map('2',1800,True),(101,'-'))

class ScoreTests(unittest.TestCase):
    def test_incomplete_score_is_never_complete_total(self):
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/'score.gz'
            with gzip.open(path,'wt') as f:f.write('# synthetic\nhm_chr\thm_pos\teffect_allele\tother_allele\teffect_weight\trsID\n1\t100\tA\tG\t0.5\trs9001\n1\t200\tC\tT\t1\trs9002\n')
            definition=[dict(id='SYNTHETIC',path=str(path),build='37',source='https://example.org')]
            partial=direct_scores(parse_genome(b'# build 37\nrs9001 1 100 AG\n'),definition,set())['scores'][0]
            self.assertEqual(partial['status'],'incomplete_coverage');self.assertIsNone(partial['weighted_sum'])
            full=direct_scores(parse_genome(b'# build 37\nrs9001 1 100 AG\nrs9002 1 200 CC\n'),definition,set())['scores'][0]
            self.assertEqual(full['weighted_sum'],2.5)
            flagged=direct_scores(parse_genome(b'# build 37\nrs9001 1 100 AG\nrs9002 1 200 CC\n'),definition,{('1',200)})['scores'][0]
            self.assertIsNone(flagged['weighted_sum'])

if __name__=='__main__':unittest.main()

class ExtraAnalysisTests(unittest.TestCase):
    def test_bed_decoding_and_empty_selection(self):
        import numpy as np
        from genetics_app.extras import bed_matrix
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/'panel'
            p.with_suffix('.fam').write_text('f a 0 0 0 -9\nf b 0 0 0 -9\nf c 0 0 0 -9\nf d 0 0 0 -9\n')
            p.with_suffix('.bim').write_text('1 rs9001 0 100 A G\n')
            p.with_suffix('.bed').write_bytes(bytes([108,27,1,0b11100100]))
            x,_,_=bed_matrix(p);self.assertEqual(x[[0,2,3],0].tolist(),[2,1,0]);self.assertTrue(np.isnan(x[1,0]))
            x,_,_=bed_matrix(p,set());self.assertEqual(x.shape,(4,0))
    def test_paternal_path_with_unassayed_ancestor(self):
        from genetics_app.extras import y_lineage
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)/'tree.json';p.write_text(json.dumps(dict(roots=[],allNodes={'1':dict(name='ROOT',parentId=None,variants=[]),'2':dict(name='TEST_BRANCH',parentId=1,variants=[dict(position=100,variant='test-marker',ancestral='A',derived='G')])})))
            data='# build 38\n'+'\n'.join(f'rs{9000+i} Y {100+i} GG' for i in range(10))
            result=y_lineage(parse_genome(data.encode()),dict(tree=str(p),source='https://example.org'))
            self.assertEqual(result['candidates'][0]['lineage'],'TEST_BRANCH');self.assertEqual(result['candidates'][0]['derived'],1)
            data=data.replace('Y 100 GG','Y 100 AA')
            self.assertEqual(y_lineage(parse_genome(data.encode()),dict(tree=str(p),source='https://example.org'))['candidates'],[])
    def test_windows_cover_targets_once(self):
        from genetics_app.setup import planned_windows
        rows=[dict(hm_chr='1',hm_pos='2000000'),dict(hm_chr='1',hm_pos='2500000'),dict(hm_chr='1',hm_pos='7000000')]
        windows=planned_windows(rows);self.assertEqual(len(windows),2)
        for r in rows:self.assertEqual(sum(w['start']<=int(r['hm_pos'])<=w['end'] for w in windows),1)

class PipelineStatusTests(unittest.TestCase):
    setUp=KnowledgeTests.setUp
    tearDown=KnowledgeTests.tearDown
    def test_missing_research_references_cannot_report_full_completion(self):
        from genetics_app.engine import analyze as basic
        from genetics_app.pipeline import full_report
        data=b'# build 37\nrs9001 1 100 AG\nrs9002 1 200 AG\n'
        profile=self.root/'pipeline.json';profile.write_text(json.dumps(dict(schema=1,knowledge=str(self.db),clinvar='unused-test-reference')))
        with patch('genetics_app.pipeline.analyze',return_value=basic(data,traits_only=True)):
            events=[]
            result=full_report(data,profile,events=events.append)
        finished=[e for e in events if e['kind']=='stage_finished']
        self.assertEqual([e['completed'] for e in finished],list(range(1,len(result['pipeline']['stages'])+1)))
        self.assertEqual(events[-1]['kind'],'complete')
        self.assertEqual(events[-1]['completed'],events[-1]['total'])
        self.assertEqual(events[-1]['status'],'completed_with_unresolved_stages')
        self.assertEqual(result['pipeline']['status'],'completed_with_unresolved_stages')
        self.assertIn('reference_not_configured',[s['status'] for s in result['pipeline']['stages']])
        self.assertEqual(result['knowledge']['counts']['genosets'],1)
    def test_keyboard_interrupt_removes_temporary_work_directory(self):
        from genetics_app.engine import analyze as basic
        from genetics_app.pipeline import full_report
        data=b'# build 37\nrs9001 1 100 AG\nrs9002 1 200 AG\n'
        profile=self.root/'cancel.json'
        profile.write_text(json.dumps(dict(schema=1,knowledge=str(self.db),clinvar='unused-test-reference',ancestry=[dict(name='Cancel fixture')])))
        work=[]
        def interrupt(genome,panel,folder,flagged):
            folder.mkdir(parents=True);(folder/'temporary-input.txt').write_text('Synthetic fixture')
            work.append(folder.parent)
            raise KeyboardInterrupt()
        with patch('genetics_app.pipeline.analyze',return_value=basic(data,traits_only=True)), patch('genetics_app.extras.ancestry',side_effect=interrupt), self.assertRaises(KeyboardInterrupt):
            full_report(data,profile)
        self.assertTrue(work);self.assertFalse(work[0].exists())

    def test_reference_change_stops_analysis(self):
        from genetics_app.pipeline import full_report
        profile=self.root/'pipeline.json';profile.write_text(json.dumps(dict(schema=1,reference_files=[dict(file=str(self.db),sha256='incorrect')])) )
        with self.assertRaisesRegex(ValueError,'checksum changed'):full_report(b'# build 37\nrs9001 1 100 AG\n',profile)
