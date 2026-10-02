"""Unit tests for the tumour-ploidy baseline and the integer SAVANA copy numbers (run: python3 -m unittest discover -s tests)."""
import math, os, sys, tempfile, unittest
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from padfoot import annot


class Rounding(unittest.TestCase):
    def test_round_half_up_is_not_bankers_rounding(self):
        self.assertEqual(annot.round_half_up(2.5), 3)
        self.assertEqual(annot.round_half_up(1.5), 2)
        self.assertEqual(annot.round_half_up(1.49), 1)
        self.assertEqual(annot.round_half_up(0.0), 0)

    def test_haplotype_baseline_is_round_of_half_the_ploidy(self):
        self.assertEqual(annot.haplotype_baseline(2.19), 1)   # SAVANA fit of P215003155
        self.assertEqual(annot.haplotype_baseline(2.29), 1)   # Wakhan fit of P215003155
        self.assertEqual(annot.haplotype_baseline(2.0), 1)
        self.assertEqual(annot.haplotype_baseline(3.0), 2)    # 1.5 rounds up
        self.assertEqual(annot.haplotype_baseline(3.8), 2)    # whole-genome doubled tumour
        self.assertEqual(annot.haplotype_baseline(5.0), 3)

    def test_cn_state(self):
        self.assertEqual(annot.cn_state(2, 1), 'AMP')
        self.assertEqual(annot.cn_state(0, 1), 'DEL')
        self.assertEqual(annot.cn_state(1, 1), 'NEUT')
        self.assertEqual(annot.cn_state(1.0, 1), 'NEUT')
        self.assertEqual(annot.cn_state(None, 1), '')
        self.assertEqual(annot.cn_state(float('nan'), 1), '')
        self.assertEqual(annot.cn_state('NA', 1), '')

    def test_round_savana_cn(self):
        self.assertEqual(annot.round_savana_cn(2.65, 0.81), (2, 1))     # OR4F5 segment: was DEL/DEL, is NEUT/NEUT at baseline 1
        self.assertEqual(annot.round_savana_cn(2.6371, 0.8015), (2, 1))
        self.assertEqual(annot.round_savana_cn(1.0, 0.0), (1, 0))       # one-copy loss, LOH
        self.assertEqual(annot.round_savana_cn(0.3, 0.1), (0, 0))       # homozygous deletion
        self.assertEqual(annot.round_savana_cn(4.0844, 1.9), (2, 2))    # balanced gain
        self.assertEqual(annot.round_savana_cn(3.0, 1.5), (2, 1))       # minor capped at floor(total/2): minor <= major
        self.assertEqual(annot.round_savana_cn(2.4, 1.3), (1, 1))

    def test_missing_minor_or_total_stays_unknown_not_zero(self):
        self.assertEqual(annot.round_savana_cn(1.1216, float('nan')), (None, None))
        self.assertEqual(annot.round_savana_cn(1.1216, None), (None, None))
        self.assertEqual(annot.round_savana_cn(float('nan'), 0.5), (None, None))

    def test_fmt_cn(self):
        self.assertEqual(annot.fmt_cn(2.0), '2')
        self.assertEqual(annot.fmt_cn(2), '2')
        self.assertEqual(annot.fmt_cn(1.84), '1.84')
        self.assertEqual(annot.fmt_cn(None), 'NA')
        self.assertEqual(annot.fmt_cn('NA'), 'NA')


class PloidyFile(unittest.TestCase):
    def _write(self, text):
        fh = tempfile.NamedTemporaryFile('w', suffix='.tsv', delete=False)
        fh.write(text)
        fh.close()
        self.addCleanup(os.remove, fh.name)
        return fh.name

    def test_savana_fitted_purity_ploidy_rank_1_row(self):
        p = self._write("purity\tploidy\tdistance\trank\n0.70\t3.95\t0.201\t2\n0.88\t2.19\t0.129\t1\n")
        self.assertAlmostEqual(annot.read_ploidy_file(p), 2.19)

    def test_wakhan_solutions_ranks_rank_1_row(self):
        p = self._write("repository_name\tdna_purity\tcell_purity\tploidy\tconfidence\tsolution_rank\n"
                        "4.58_0.93_0.80\t0.94\t0.93\t4.58\t0.80\t2\n2.29_0.93_0.91\t0.94\t0.93\t2.29\t0.91\t1\n")
        self.assertAlmostEqual(annot.read_ploidy_file(p), 2.29)

    def test_plain_table_without_rank_uses_first_row(self):
        p = self._write("ploidy\n3.4\n2.0\n")
        self.assertAlmostEqual(annot.read_ploidy_file(p), 3.4)

    def test_missing_column_or_bad_value_raises(self):
        with self.assertRaises(ValueError):
            annot.read_ploidy_file(self._write("purity\tdistance\n0.9\t0.1\n"))
        with self.assertRaises(ValueError):
            annot.read_ploidy_file(self._write("ploidy\trank\n0\t1\n"))
        with self.assertRaises(ValueError):
            annot.read_ploidy_file(self._write("ploidy\trank\nNA\t1\n"))


class MeanAndResolve(unittest.TestCase):
    def test_mean_cn_is_length_weighted_and_skips_unknown(self):
        segs = [(2.0, 1, 100), (4.0, 101, 400), (None, 401, 1000), (float('nan'), 1001, 2000)]
        self.assertAlmostEqual(annot.mean_cn(segs), (2.0 * 100 + 4.0 * 300) / 400)
        self.assertTrue(math.isnan(annot.mean_cn([(None, 1, 10)])))

    def test_resolve_ploidy_prefers_the_given_value(self):
        ploidy, base = annot.resolve_ploidy(3.8, 2.0, 'test')
        self.assertEqual((ploidy, base), (3.8, 2))

    def test_resolve_ploidy_falls_back_to_the_estimate_with_a_warning(self):
        with self.assertLogs(annot.logger, level='WARNING') as cm:
            ploidy, base = annot.resolve_ploidy(None, 2.4, 'SAVANA')
        self.assertEqual((ploidy, base), (2.4, 1))
        self.assertTrue(any('No tumour ploidy supplied' in m for m in cm.output))
        ploidy, base = annot.resolve_ploidy(None, float('nan'), 'x')
        self.assertEqual((ploidy, base), (2.0, 1))


SAVANA_HEADER = "chromosome\tstart\tend\tsegment_id\tbin_count\tsum_of_bin_lengths\tweight\tcopyNumber\tminorAlleleCopyNumber\tmeanBAF\tno_hetSNPs\n"


class GetSavanaCNA(unittest.TestCase):
    def _table(self, rows):
        fh = tempfile.NamedTemporaryFile('w', suffix='.tsv', delete=False)
        fh.write(SAVANA_HEADER)
        for i, (chrom, s, e, cn, minor) in enumerate(rows):
            fh.write(f"{chrom}\t{s}\t{e}\t{chrom}_seg{i}\t10\t{e - s + 1}\t10.0\t{cn}\t{minor}\t0.5\t100\n")
        fh.close()
        self.addCleanup(os.remove, fh.name)
        return fh.name

    def test_integer_alleles_and_labels_against_the_fitted_ploidy(self):
        p = self._table([('chr19', 1, 1000000, 2.65, 0.81),          # diploid with noise -> 2/1 -> NEUT/NEUT
                         ('chr19', 1000001, 2000000, 1.05, 0.02),    # one-copy loss -> 1/0 -> NEUT/DEL, LOH
                         ('chr19', 2000001, 3000000, 3.1, 0.95),     # gain of one allele -> 2/1 -> AMP/NEUT
                         ('chr19', 3000001, 4000000, 4.0844, '')])   # no het SNPs -> unknown/unknown
        cnas, baselines = annot.get_savana_CNA(p, [], ploidy=2.19)
        self.assertEqual(baselines, [1, 1])
        hp1 = {(c.pos_1): (c.cn, c.dir1, c.LOH) for c in cnas[('chr19', 1)]}
        hp2 = {(c.pos_1): (c.cn, c.dir1, c.LOH) for c in cnas[('chr19', 2)]}
        self.assertEqual(hp1[1], (2, 'AMP', False))
        self.assertEqual(hp2[1], (1, 'NEUT', False))
        self.assertEqual((hp1[1000001], hp2[1000001]), ((1, 'NEUT', True), (0, 'DEL', True)))
        self.assertEqual((hp1[2000001][:2], hp2[2000001][:2]), ((2, 'AMP'), (1, 'NEUT')))
        self.assertEqual((hp1[3000001], hp2[3000001]), ((None, '', False), (None, '', False)))

    def test_fallback_estimate_is_the_mean_total_copy_number(self):
        p = self._table([('chr1', 1, 1000, 2.0, 1.0), ('chr1', 1001, 4000, 4.0, 2.0)])
        with self.assertLogs(annot.logger, level='WARNING'):
            cnas, baselines = annot.get_savana_CNA(p, [], ploidy=None)   # mean total = 3.5 -> baseline 2
        self.assertEqual(baselines, [2, 2])
        self.assertEqual([c.dir1 for c in cnas[('chr1', 1)]], ['DEL', 'NEUT'])

    def test_first_segment_of_the_real_table_is_not_a_double_deletion(self):
        # P215003155 chr1:10001-770000: copyNumber 2.6371, minorAlleleCopyNumber 0.8015, SAVANA ploidy 2.19
        p = self._table([('chr1', 10001, 770000, 2.6371, 0.8015)])
        cnas, _ = annot.get_savana_CNA(p, [], ploidy=2.19)
        self.assertEqual((cnas[('chr1', 1)][0].dir1, cnas[('chr1', 2)][0].dir1), ('AMP', 'NEUT'))


WAKHAN_HEADER = """##fileformat=VCFv4.2
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="x">
##INFO=<ID=END,Number=1,Type=Integer,Description="x">
##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="x">
##INFO=<ID=BPS,Number=0,Type=String,Description="x">
##FORMAT=<ID=GT,Number=1,Type=String,Description="x">
##FORMAT=<ID=TCN,Number=1,Type=Float,Description="x">
##FORMAT=<ID=CN1,Number=1,Type=Float,Description="x">
##FORMAT=<ID=CN2,Number=1,Type=Float,Description="x">
##FORMAT=<ID=CNQ1,Number=1,Type=Float,Description="x">
##FORMAT=<ID=CNQ2,Number=1,Type=Float,Description="x">
##FORMAT=<ID=COV1,Number=1,Type=Float,Description="x">
##FORMAT=<ID=COV2,Number=1,Type=Float,Description="x">
##contig=<ID=chr1,length=10000>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample
"""


class GetCNAWakhan(unittest.TestCase):
    def _vcf(self, records):
        fh = tempfile.NamedTemporaryFile('w', suffix='.vcf', delete=False)
        fh.write(WAKHAN_HEADER)
        for chrom, s, e, cn1, cn2 in records:
            fh.write(f"{chrom}\t{s}\twakhan:x:{chrom}:{s}-{e}\tN\t<DUP>\t1000\tPASS\tSVTYPE=CNV;SVLEN={e - s};END={e};BPS=\t"
                     f"GT:TCN:CN1:CN2:CNQ1:CNQ2:COV1:COV2\t0/1:{cn1 + cn2}:{cn1}:{cn2}:0.9:0.9:{15.0 * cn1}:{15.0 * cn2}\n")
        fh.close()
        self.addCleanup(os.remove, fh.name)
        return fh.name

    def test_diploid_tumour_fillers_are_neutral_and_gains_amplified(self):
        p = self._vcf([('chr1', 2001, 4000, 2.0, 1.0)])
        cnas, baselines = annot.get_CNA(p, [], ploidy=2.29)
        self.assertEqual(baselines, [1, 1])
        self.assertEqual([(c.cn, c.dir1) for c in cnas[('chr1', 1)]], [(1.0, 'NEUT'), (2.0, 'AMP'), (1.0, 'NEUT')])
        self.assertEqual([(c.cn, c.dir1) for c in cnas[('chr1', 2)]], [(1.0, 'NEUT'), (1.0, 'NEUT'), (1.0, 'NEUT')])

    def test_whole_genome_doubled_tumour_uses_the_fitted_ploidy(self):
        # ploidy 3.8 -> baseline 2: the CN 2 segment is neutral, the unlisted (reference 1/1) padding is a relative loss
        p = self._vcf([('chr1', 2001, 4000, 2.0, 2.0)])
        cnas, baselines = annot.get_CNA(p, [], ploidy=3.8)
        self.assertEqual(baselines, [2, 2])
        self.assertEqual([c.dir1 for c in cnas[('chr1', 1)]], ['DEL', 'NEUT', 'DEL'])

    def test_fallback_estimate_sums_both_haplotype_means(self):
        p = self._vcf([('chr1', 1, 10000, 2.0, 1.0)])   # whole contig listed as 2/1 -> mean total 3.0 -> baseline 2
        with self.assertLogs(annot.logger, level='WARNING'):
            cnas, baselines = annot.get_CNA(p, [], ploidy=None)
        self.assertEqual(baselines, [2, 2])
        self.assertEqual([c.dir1 for c in cnas[('chr1', 1)]], ['NEUT'])
        self.assertEqual([c.dir1 for c in cnas[('chr1', 2)]], ['DEL'])


class AnnotCNAsUnknown(unittest.TestCase):
    def test_unknown_copy_number_gives_NA_and_no_label(self):
        from collections import defaultdict
        cnas = {('chr1', 1): [annot.CNA('chr1', 1, 5000, None, 1, False, '', ''), annot.CNA('chr1', 5001, 10000, 2, 1, False, '', '')],
                ('chr1', 2): [annot.CNA('chr1', 1, 5000, None, 2, False, '', ''), annot.CNA('chr1', 5001, 10000, 1, 2, False, '', '')]}
        for key in cnas:
            for c in cnas[key]:
                c.dir1 = annot.cn_state(c.cn, 1)
        genes = {'chr1': [['GENE_A', 'GENE_B'], [1000, 6000], [2000, 7000]]}
        by_gene = defaultdict(list)
        annot.annot_CNAs(genes, cnas, [1, 1], by_gene)
        self.assertEqual((by_gene['GENE_A'].CN, by_gene['GENE_A'].CN_impact), (['NA', 'NA'], ['', '']))
        self.assertEqual((by_gene['GENE_B'].CN, by_gene['GENE_B'].CN_impact), ([2, 1], ['AMP', 'NEUT']))
        row = by_gene['GENE_A'].to_str().split('\t')
        self.assertEqual((row[4], row[5], row[8], row[9]), ('', 'NA', '', 'NA'))


if __name__ == '__main__':
    unittest.main()
