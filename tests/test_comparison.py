import unittest
from tools.head_to_head import compare_rows


class ComparisonTests(unittest.TestCase):
    def report(self):
        return {'knowledge':{'findings':[dict(id='rs1',kind='variant',source_genotype='(A;G)',raw_ids=['rs1','i2'],summary='Example',has_genotype_summary=True,repute='Good',magnitude=2)],'unresolved':[],'genoset_audit':[],'counts':{}}}

    def test_coordinate_colocation_is_not_alias_evidence(self):
        rows=[dict(rsnum='i2',geno='(A;G)',genosummary='Example')]
        self.assertEqual(compare_rows(rows,self.report(),{})['counts'],{'unresolved':1})
        self.assertEqual(compare_rows(rows,self.report(),{'i2':'rs1'})['counts'],{'explicit_alias_represented':1})

    def test_discordance_cannot_be_hidden_by_summary_match(self):
        rows=[dict(rsnum='rs1',geno='(A;A)',genosummary='Example')]
        r=compare_rows(rows,self.report(),{})
        self.assertEqual(r['counts'],{'genotype_disagreement':1})
        self.assertEqual(r['content_agreement'],{})

    def test_allele_order_and_markup_normalized_but_real_text_difference_retained(self):
        rows=[dict(rsnum='rs1',geno='(G;A)',genosummary='<b>Example</b>',magnitude=2,repute='Good')]
        self.assertEqual(compare_rows(rows,self.report(),{})['content_agreement']['baseline_explanations_same'],1)
        rows[0]['genosummary']='Different interpretation'
        self.assertEqual(compare_rows(rows,self.report(),{})['content_agreement']['baseline_explanations_changed'],1)

    def test_sequence_notation_requires_explicit_mapping_and_matching_copies(self):
        r=self.report();f=r['knowledge']['findings'][0];f['source_genotype']='(-;CT)'
        rows=[dict(rsnum='rs1',geno='(D;I)')]
        self.assertEqual(compare_rows(rows,r,{})['counts'],{'genotype_disagreement':1})
        f['sequence_resolution']={'definition':{'source_insertion':'CT'}}
        self.assertEqual(compare_rows(rows,r,{})['counts'],{'public_sequence_notation_equivalent':1})
        rows[0]['geno']='(I;I)'
        self.assertEqual(compare_rows(rows,r,{})['counts'],{'genotype_disagreement':1})

    def test_strand_notation_is_separate_and_needs_a_reviewed_definition(self):
        r=self.report();f=r['knowledge']['findings'][0];f.update(raw_genotype='CC',source_genotype='(G;G)')
        rows=[dict(rsnum='rs1',geno='(C;C)')]
        self.assertEqual(compare_rows(rows,r,{})['counts'],{'genotype_disagreement':1})
        f['curated_sequence_evidence']={'orientation':'minus'}
        self.assertEqual(compare_rows(rows,r,{})['counts'],{'public_strand_notation_equivalent':1})
        rows[0]['geno']='(T;T)'
        self.assertEqual(compare_rows(rows,r,{})['counts'],{'genotype_disagreement':1})
        f['curated_sequence_evidence']['orientation']='plus';rows[0]['geno']='(C;C)'
        self.assertEqual(compare_rows(rows,r,{})['counts'],{'genotype_disagreement':1})
