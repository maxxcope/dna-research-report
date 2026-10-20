"""Synthetic checks of the audit prototype, not clinical validation."""
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'tools'))
from snpedia_rules import UnsupportedRule, evaluate, evaluate_all, parse_rule
from audit_snpedia_parity import marker_metadata, oriented_calls
from genetics_app.parser import parse_genome


class RuleTests(unittest.TestCase):
    def run_rule(self, source, calls):
        return evaluate(parse_rule(source), calls, lambda _: None)

    def test_missing_under_not_does_not_become_normal(self):
        self.assertIsNone(self.run_rule('not(rs901(T;T))', {}))
        self.assertIs(self.run_rule('not(rs901(T;T))', {'rs901': ('A', 'A')}), True)

    def test_threshold_distinguishes_unknown_from_false(self):
        rule = 'atleast(2,rs901(A),rs902(A),rs903(A))'
        self.assertIsNone(self.run_rule(rule, {'rs901': ('A', 'G')}))
        self.assertIs(self.run_rule(rule, {'rs901': ('G', 'G'), 'rs902': ('G', 'G')}), False)
        self.assertIs(self.run_rule(rule, {'rs901': ('A', 'G'), 'rs902': ('A', 'A')}), True)

    def test_all_operator_truth_tables(self):
        values = {True: ('A', 'A'), False: ('G', 'G'), None: None}
        for left in values:
            for right in values:
                calls = {key: values[value] for key, value in [('rs901', left), ('rs902', right)] if value is not None}
                expected_and = False if left is False or right is False else None if left is None or right is None else True
                expected_or = True if left is True or right is True else None if left is None or right is None else False
                self.assertIs(self.run_rule('and(rs901(A),rs902(A))', calls), expected_and)
                self.assertIs(self.run_rule('or(rs901(A),rs902(A))', calls), expected_or)
                expected_not = None if expected_or is None else not expected_or
                self.assertIs(self.run_rule('not(rs901(A),rs902(A))', calls), expected_not)

    def test_allele_order_and_single_allele(self):
        self.assertTrue(self.run_rule('rs901(G;A)', {'rs901': ('A', 'G')}))
        self.assertFalse(self.run_rule('rs901(G;G)', {'rs901': ('A', 'G')}))
        self.assertTrue(self.run_rule('rs901(G)', {'rs901': ('A', 'G')}))

    def test_real_grammar_comments_and_trailing_comma(self):
        rule = ' # description\n <!-- hidden --> and(rs901(A;G), or(rs902(T),rs903(T),),)'
        self.assertTrue(self.run_rule(rule, {'rs901': ('G', 'A'), 'rs903': ('T', 'T')}))

    def test_forward_references_and_cycles(self):
        outcomes, errors = evaluate_all({'gs901': 'gs902', 'gs902': 'rs901(A)'}, {'rs901': ('A', 'A')})
        self.assertEqual(outcomes, {'gs902': True, 'gs901': True})
        self.assertFalse(errors)
        outcomes, _ = evaluate_all({'gs901': 'gs902', 'gs902': 'gs901'}, {})
        self.assertTrue(all(value is None for value in outcomes.values()))
        outcomes, _ = evaluate_all({'gs901': 'or(gs901,rs901(A))'}, {'rs901': ('A', 'A')})
        self.assertIsNone(outcomes['gs901'])

    def test_untrusted_or_malformed_rules_never_execute(self):
        for source in ["__import__('os').system('anything')", 'rs901(A); rs902(T)', 'and()',
                       'atleast(-1,rs901(A))', 'atleast(2,rs901(A))', 'rs901(N)', 'not()']:
            with self.subTest(source=source), self.assertRaises(UnsupportedRule):
                parse_rule(source)


class SourceMatchingTests(unittest.TestCase):
    def test_metadata_not_first_template_and_inline_fields(self):
        source = '{{PMID|123}}\n{{Rsnum|rsid=901|Orientation=plus|StabilizedOrientation=minus}}'
        self.assertEqual(marker_metadata(source)['stabilizedorientation'], 'minus')

    def test_strand_and_exclusions(self):
        genome = parse_genome(b'# build 37\nrs901 1 100 AG\nrs902 1 200 AA\nrs903 1 300 AA\nrs904 1 400 II\n')
        snps = {'rs901': '{{Rsnum|StabilizedOrientation=minus}}',
                'rs902': '{{Rsnum|Orientation=plus|Flip37=true}}',
                'rs903': '{{Rsnum|rsid=903}}',
                'rs904': '{{Rsnum|Orientation=plus}}'}
        calls, reasons, _ = oriented_calls(genome, snps)
        self.assertEqual(calls, {'rs901': ('C', 'T')})
        self.assertEqual(reasons['rs902'], 'build_flip_unresolved')
        self.assertEqual(reasons['rs903'], 'orientation_unknown')
        self.assertEqual(reasons['rs904'], 'indel_code_unresolved')

    def test_alias_conflict_withheld_and_direct_not_overwritten(self):
        genome = parse_genome(b'# build 37\ni901 1 100 AG\ni902 1 200 AA\nrs903 1 300 GG\ni903 1 400 AA\n')
        snps = {'i901': '{{23andMe SNP|rsid=901}}', 'i902': '{{23andMe SNP|rsid=901}}',
                'rs901': '{{Rsnum|Orientation=plus}}', 'rs903': '{{Rsnum|Orientation=plus}}',
                'i903': '{{23andMe SNP|rsid=903}}'}
        calls, reasons, _ = oriented_calls(genome, snps, True)
        self.assertNotIn('rs901', calls)
        self.assertEqual(reasons['rs901'], 'alias_probe_unresolved')
        self.assertEqual(calls['rs903'], ('G', 'G'))


if __name__ == '__main__':
    unittest.main()
