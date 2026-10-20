import csv
import io
import json
from pathlib import Path
import tempfile
import threading
import unittest
from unittest.mock import patch
import urllib.request
import urllib.error
from genetics_app.parser import parse_genome, InputError
from genetics_app.engine import analyze, classification_category
from genetics_app.references import import_clinvar, import_gwas, metadata
from genetics_app.report import report_lines


def raw(rows,build='37'):
    return (f'# build {build}\n'+rows).encode()

class ParserTests(unittest.TestCase):
    def test_four_column_and_bom(self):
        g=parse_genome(b'\xef\xbb\xbf'+raw('rs1 1 100 GA\nrs2 MT 20 T\n'))
        self.assertEqual(g.locus('1',100)[:2],('called','AG'))
        self.assertEqual(g.locus('MT',20)[:2],('called','T'))
    def test_five_column_missing(self):
        g=parse_genome(raw('rsid chromosome position allele1 allele2\nrs1 1 100 A G\nrs2 2 200 0 0\n'))
        self.assertEqual(g.locus('2',200)[0],'no_call')
    def test_csv(self):
        self.assertEqual(parse_genome(raw('"RSID","CHROMOSOME","POSITION","RESULT"\n"rs1","1","100","AG"\n')).qc['rows'],1)
    def test_build_required(self):
        with self.assertRaises(InputError):parse_genome(b'rs1 1 100 AA')
        self.assertEqual(parse_genome(b'rs1 1 100 AA',build='38').build,'38')
    def test_build_conflict(self):
        with self.assertRaises(InputError):parse_genome(raw('rs1 1 100 AA'),'38')
    def test_unsupported_build(self):
        with self.assertRaises(InputError):parse_genome(raw('rs1 1 100 AA','36'),'37')
    def test_ambiguous_header(self):
        with self.assertRaises(InputError):parse_genome(b'# build 37 and GRCh38\nrs1 1 100 AA')
    def test_duplicates_fail_closed(self):
        g=parse_genome(raw('rs1 1 100 AA\nrs1 1 101 AA\nrs2 2 200 AA\nrs3 2 200 CC\n'))
        self.assertEqual(g.locus('1',100)[0],'conflicting_marker_id')
        self.assertEqual(g.locus('2',200)[0],'conflicting_probes')
    def test_identical_duplicate(self):
        g=parse_genome(raw('rs1 1 100 AG\nrs1 1 100 GA\n'))
        self.assertEqual(g.locus('1',100)[0],'called')
    def test_no_call_and_called_probe(self):
        g=parse_genome(raw('rs1 1 100 --\nrs2 1 100 AG\n'))
        self.assertEqual(g.locus('1',100)[1],'AG')
    def test_indel_and_autosomal_haploid_withheld(self):
        g=parse_genome(raw('rs1 1 100 II\nrs2 2 200 A\n'))
        self.assertEqual(g.locus('1',100)[0],'indel_code_unresolved')
        self.assertEqual(g.locus('2',200)[0],'unexpected_autosomal_ploidy')
    def test_bad_rows(self):
        for row in ['rs1 25 100 AA','rs1 1 0 AA','rs1 1 100 A-','rs1 1 100 AG X X','<script> 1 100 AA']:
            with self.subTest(row=row),self.assertRaises(InputError):parse_genome(raw(row))
    def test_synthetic_label(self):
        self.assertTrue(parse_genome(b'# Synthetic example\n'+raw('rs1 1 100 AA')).qc['synthetic_input_declared'])
    def test_forward_only(self):
        with self.assertRaises(InputError):parse_genome(raw('rs1 1 100 AA'),strand='reverse')

class EngineTests(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.root=Path(self.tmp.name)
        self.vcf=self.root/'ref.vcf';self.db=self.root/'clinvar37.sqlite'
        self.vcf.write_text('''##fileformat=VCFv4.2
##fileDate=2026-10-01
##reference=GRCh37
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
1\t100\t1\tA\tG\t.\t.\tRS=1;CLNSIG=Pathogenic;CLNREVSTAT=reviewed_by_expert_panel;GENEINFO=TEST:1;CLNDN=Mock_condition;ORIGIN=1
1\t101\t2\tA\tG\t.\t.\tRS=2;CLNSIG=Benign;CLNREVSTAT=criteria_provided,_single_submitter;ORIGIN=1
1\t102\t3\tA\tG\t.\t.\tRS=3;CLNSIG=Conflicting_classifications_of_pathogenicity;CLNREVSTAT=criteria_provided,_conflicting_classifications;ORIGIN=1
1\t103\t4\tA\tG\t.\t.\tRS=4;CLNSIG=Uncertain_significance;ORIGIN=1
1\t104\t5\tA\tAT\t.\t.\tRS=5;CLNSIG=Pathogenic;ORIGIN=1
1\t105\t6\tA\tG,T\t.\t.\tRS=6;CLNSIG=Pathogenic;ORIGIN=1
1\t106\t7\tA\tG\t.\t.\tRS=7;CLNSIG=Pathogenic;ORIGIN=2
1\t107\t8\tA\tG\t.\t.\tRS=8;CLNSIG=drug_response;ORIGIN=1
''')
        import_clinvar(self.vcf,self.db,'37')
    def tearDown(self):self.tmp.cleanup()
    def test_low_penetrance_and_combined_classes(self):
        self.assertEqual(classification_category('Pathogenic/Likely_pathogenic'),'Review with a professional')
        self.assertEqual(classification_category('Pathogenic,_low_penetrance'),'Review with a professional')
        self.assertEqual(classification_category('Benign/Likely_benign'),'Benign / likely benign')
    def test_matching_and_classification(self):
        r=analyze(raw('rs1 1 100 AG\nrs2 1 101 AA\nrs3 1 102 GG\nrs4 1 103 --\nrs5 1 104 AA\nrs6 1 105 AG\nrs7 1 106 AG\nrs8 1 107 TT\n'),clinvar=self.db)
        self.assertEqual(len(r['clinvar']),2)
        self.assertEqual(r['clinvar'][0]['category'],'Review with a professional')
        self.assertEqual(r['clinvar'][0]['stars'],3)
        self.assertEqual(r['clinvar'][1]['category'],'Conflicting evidence')
        counts=r['clinvar_counts']
        for k,v in {'reference_only':1,'no_call':1,'non_biallelic_snv_unresolved':2,'somatic_only_unresolved':1,'allele_mismatch':1}.items():self.assertEqual(counts[k],v)
    def test_no_complement_guess(self):
        r=analyze(raw('rs1 1 100 CT\n'),clinvar=self.db)
        self.assertEqual(r['clinvar'],[])
        self.assertEqual(r['exclusions'][0]['reason'],'allele_mismatch')
    def test_ref_build_mismatch(self):
        with self.assertRaises(ValueError):analyze(raw('rs1 1 100 AG','38'),clinvar=self.db)
    def test_import_build_mismatch_atomic(self):
        before=self.db.read_bytes()
        with self.assertRaises(ValueError):import_clinvar(self.vcf,self.db,'38')
        self.assertEqual(before,self.db.read_bytes())
    def test_unsupported_data_never_normal(self):
        r=analyze(raw('rs1 1 100 II\nrs2 1 101 --\n'),clinvar=self.db)
        self.assertEqual(r['clinvar'],[])
        self.assertEqual(len(r['exclusions']),2)
    def test_traits_expected_outputs_and_forward_alleles(self):
        r=analyze(raw('rs4988235 2 136608646 AG\nrs671 12 112241766 GG\nrs17822931 16 48258198 TT\n'),traits_only=True)
        self.assertIn('One copy',r['traits'][0]['explanation'])
        self.assertIn('not present',r['traits'][1]['explanation'])
        self.assertIn('dry earwax',r['traits'][2]['explanation'])
    def test_traits_absent_wrong_coord_alias(self):
        r=analyze(raw('rs4988235 2 1 AA\ni123 16 48258198 TT\n'),traits_only=True)
        self.assertEqual(r['traits'][0]['status'],'not_assayed')
        self.assertEqual(r['traits'][2]['status'],'marker_identity_unverified')
    def test_traits_grch38(self):
        r=analyze(raw('rs17822931 16 48224287 CT','38'),traits_only=True)
        self.assertIn('wet earwax',r['traits'][2]['explanation'])
    def test_reference_required_by_default(self):
        with self.assertRaises(ValueError):analyze(raw('rs1 1 100 AA'))
    def test_offline_analysis(self):
        with patch('socket.socket',side_effect=AssertionError('Network forbidden')):
            r=analyze(raw('rs1 1 100 AG'),clinvar=self.db)
        self.assertEqual(len(r['clinvar']),1)
    def test_gwas_leads_not_risk_predictions(self):
        p=self.root/'gwas.tsv';fields=['SNPS','P-VALUE','PUBMEDID','DISEASE/TRAIT','STRONGEST SNP-RISK ALLELE']
        with p.open('w') as f:
            w=csv.DictWriter(f,fieldnames=fields,delimiter='\t');w.writeheader();w.writerows([
                dict(zip(fields,['rs1','1e-9','123','Mock trait','rs1-G'])),dict(zip(fields,['rs2','0.1','124','Rejected','rs2-G']))])
        gdb=self.root/'gwas.sqlite';import_gwas(p,gdb)
        r=analyze(raw('rs1 1 100 AA'),clinvar=self.db,gwas=gdb)
        self.assertEqual(len(r['gwas']),1)
        self.assertIn('no personal risk direction',r['gwas'][0]['interpretation'])
        self.assertEqual(r['clinvar'],[])

if __name__=='__main__':unittest.main()
