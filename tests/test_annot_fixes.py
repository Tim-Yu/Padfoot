"""Unit tests for pad_unlisted_segments() and get_repeat() (run: python3 -m unittest discover -s tests)."""
import os, sys, tempfile, unittest
from collections import defaultdict
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from padfoot import annot


class PadUnlistedSegments(unittest.TestCase):
    def test_ends_gap_and_untouched_chromosome(self):
        hpls = {'chr1': [[0.0, 2.0], [100, 301], [200, 400]]}
        out = annot.pad_unlisted_segments({'chr1': 1000, 'chr2': 500}, hpls)
        self.assertEqual(out['chr1'], [[1.0, 0.0, 1.0, 2.0, 1.0], [1, 100, 201, 301, 401], [99, 200, 300, 400, 1000]])
        self.assertEqual(out['chr2'], [[1.0], [1], [500]])

    def test_adjacent_segments_with_different_cn_are_not_merged(self):
        hpls = {'chr1': [[2.0, 0.0], [1, 201], [200, 1000]]}
        out = annot.pad_unlisted_segments({'chr1': 1000}, hpls)
        self.assertEqual(out['chr1'], [[2.0, 0.0], [1, 201], [200, 1000]])

    def test_unsorted_input_and_passthrough_of_contigs_without_length(self):
        hpls = {'chr1': [[2.0, 0.0], [500, 10], [600, 20]], 'chrUn': [[3.0], [1], [5]]}
        out = annot.pad_unlisted_segments({'chr1': 700}, hpls)
        self.assertEqual(out['chr1'][1], [1, 10, 21, 500, 601])
        self.assertEqual(out['chr1'][0], [1.0, 0.0, 1.0, 2.0, 1.0])
        self.assertEqual(out['chrUn'], [[3.0], [1], [5]])

    def test_lookup_no_longer_wraps_to_the_last_segment(self):
        # the failure mode: a gene before the first listed segment used to hit index -1
        import bisect
        hpls = {'chr1': [[0.0], [120000001], [151000000]]}
        out = annot.pad_unlisted_segments({'chr1': 248956422}, hpls)
        cn, starts, _ = out['chr1']
        self.assertEqual(cn[bisect.bisect_right(starts, 55419) - 1], 1.0)        # OR4F5-like gene -> neutral
        self.assertEqual(cn[bisect.bisect_right(starts, 130000000) - 1], 0.0)    # inside the loss
        self.assertEqual(cn[bisect.bisect_right(starts, 200000000) - 1], 1.0)    # after the loss -> neutral


WAKHAN_HEADER = """##fileformat=VCFv4.2
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="x">
##INFO=<ID=END,Number=1,Type=Integer,Description="x">
##INFO=<ID=BPS,Number=0,Type=String,Description="x">
##FORMAT=<ID=GT,Number=1,Type=String,Description="x">
##FORMAT=<ID=TCN,Number=1,Type=Float,Description="x">
##FORMAT=<ID=CN1,Number=1,Type=Float,Description="x">
##FORMAT=<ID=CN2,Number=1,Type=Float,Description="x">
##FORMAT=<ID=CNQ1,Number=1,Type=Float,Description="x">
##FORMAT=<ID=CNQ2,Number=1,Type=Float,Description="x">
##FORMAT=<ID=COV1,Number=1,Type=Float,Description="x">
##FORMAT=<ID=COV2,Number=1,Type=Float,Description="x">
##contig=<ID=chr3,length=198295559>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample
"""


def wakhan_record(chrom, pos, end, cn1, cn2):
    return (f"{chrom}\t{pos}\twakhan:X:{chrom}:{pos}-{end}\tN\t<DEL>\t1000\tPASS\tSVTYPE=CNV;END={end};BPS=\t"
            f"GT:TCN:CN1:CN2:CNQ1:CNQ2:COV1:COV2\t1/1:{cn1 + cn2}:{cn1}:{cn2}:0.9:0.9:{15 * cn1}:{15 * cn2}\n")


class GetCNA(unittest.TestCase):
    def test_same_cn_records_separated_by_a_gap_stay_separate_and_the_gap_is_neutral(self):
        # the real chr3 case: two 0/0 losses 69 Mb apart were merged into one 87-162.9 Mb loss
        with tempfile.TemporaryDirectory() as d:
            vcf = os.path.join(d, 'w.vcf')
            with open(vcf, 'w') as fh:
                fh.write(WAKHAN_HEADER)
                fh.write(wakhan_record('chr3', 44699740, 44700555, 2.0, 1.0))
                fh.write(wakhan_record('chr3', 87000001, 94000000, 0.0, 0.0))
                fh.write(wakhan_record('chr3', 162807955, 162908546, 0.0, 0.0))
            cnas, ploidy = annot.get_CNA(vcf, [])
        segs = [(c.pos_1, c.pos_2, c.cn) for c in cnas[('chr3', 1)]]
        self.assertEqual(segs, [(1, 44699739, 1.0), (44699740, 44700555, 2.0), (44700556, 87000000, 1.0),
                                (87000001, 94000000, 0.0), (94000001, 162807954, 1.0),
                                (162807955, 162908546, 0.0), (162908547, 198295559, 1.0)])
        self.assertEqual(ploidy, [1, 1])

    def test_records_are_sorted_before_segments_are_built(self):
        with tempfile.TemporaryDirectory() as d:
            vcf = os.path.join(d, 'w.vcf')
            with open(vcf, 'w') as fh:
                fh.write(WAKHAN_HEADER)
                fh.write(wakhan_record('chr3', 162807955, 162908546, 0.0, 0.0))   # out of order on purpose
                fh.write(wakhan_record('chr3', 87000001, 94000000, 0.0, 0.0))
                fh.write(wakhan_record('chr3', 44699740, 44700555, 2.0, 1.0))
            cnas, _ = annot.get_CNA(vcf, [])
        self.assertEqual([(c.pos_1, c.pos_2, c.cn) for c in cnas[('chr3', 1)]],
                         [(1, 44699739, 1.0), (44699740, 44700555, 2.0), (44700556, 87000000, 1.0),
                          (87000001, 94000000, 0.0), (94000001, 162807954, 1.0),
                          (162807955, 162908546, 0.0), (162908547, 198295559, 1.0)])

    def test_adjacent_same_cn_records_are_still_merged(self):
        with tempfile.TemporaryDirectory() as d:
            vcf = os.path.join(d, 'w.vcf')
            with open(vcf, 'w') as fh:
                fh.write(WAKHAN_HEADER)
                fh.write(wakhan_record('chr3', 1000001, 2000000, 0.0, 1.0))
                fh.write(wakhan_record('chr3', 2000001, 3000000, 0.0, 1.0))
                fh.write(wakhan_record('chr3', 5000001, 6000000, 2.0, 1.0))   # get_CNA needs one segment with CN1 > 0
            cnas, _ = annot.get_CNA(vcf, [])
        self.assertEqual([(c.pos_1, c.pos_2, c.cn) for c in cnas[('chr3', 1)]],
                         [(1, 1000000, 1.0), (1000001, 3000000, 0.0), (3000001, 5000000, 1.0), (5000001, 6000000, 2.0),
                          (6000001, 198295559, 1.0)])


class SV:  # minimal stand-in for annot.SV
    def __init__(self):
        self.repeat = []


HEADER = "   SW   perc perc perc  query      position in query    matching repeat        position in repeat\nscore   div. del. ins.  sequence   begin end   (left)   repeat   class/family begin  end    (left)  ID\n\n"


class GetRepeat(unittest.TestCase):
    def run_in_tmp(self, out_text, ids):
        cwd = os.getcwd()
        with tempfile.TemporaryDirectory() as d:
            os.chdir(d)
            try:
                with open('temp_ins.fa.out', 'w') as fh:
                    fh.write(out_text)
                svls = {i: SV() for i in ids}
                annot.get_repeat(svls)
                return {i: svls[i].repeat for i in ids}
            finally:
                os.chdir(cwd)

    def test_last_query_is_summarised(self):
        text = HEADER + \
            " 2000 10.0 0.0 0.0 INS_A 1 300 (0) + AluSx SINE/Alu 2 301 (11) 1\n" + \
            " 1500 12.0 0.0 0.0 INS_B 1 280 (20) + L1PA2 LINE/L1 1 280 (5000) 2\n"
        rep = self.run_in_tmp(text, ['INS_A', 'INS_B'])
        self.assertEqual(rep['INS_A'], [('SINE/Alu', 299)])
        self.assertEqual(rep['INS_B'], [('LINE/L1', 279)])   # was silently dropped before the fix

    def test_hits_of_one_query_split_over_two_blocks_are_summed(self):
        # RepeatMasker lists INS_A's hits in two non-adjacent blocks; each alone is below 80 %, together 530/545
        text = HEADER + \
            " 108 0.0 0.0 0.0 INS_A 419 535 (10) + (T)n Simple_repeat 1 117 (0) 30\n" + \
            " 2000 10.0 0.0 0.0 INS_B 1 300 (0) + AluSx SINE/Alu 2 301 (11) 1\n" + \
            " 458 0.5 0.7 0.0 INS_A 9 418 (127) + (ATTCC)n Simple_repeat 1 413 (0) 46\n"
        rep = self.run_in_tmp(text, ['INS_A', 'INS_B'])
        self.assertEqual(rep['INS_A'], [('Simple_repeat', 525)])
        self.assertEqual(rep['INS_B'], [('SINE/Alu', 299)])

    def test_single_query_and_80_percent_rule(self):
        text = HEADER + " 2000 10.0 0.0 0.0 INS_A 1 300 (0) + AluSx SINE/Alu 2 301 (11) 1\n"
        self.assertEqual(self.run_in_tmp(text, ['INS_A'])['INS_A'], [('SINE/Alu', 299)])
        text = HEADER + " 300 10.0 0.0 0.0 INS_A 1 100 (200) + AluSx SINE/Alu 2 101 (211) 1\n"
        self.assertEqual(self.run_in_tmp(text, ['INS_A'])['INS_A'], [])   # 100/300 < 80 %

    def test_no_hit_notice_and_unknown_id_do_not_crash(self):
        text = "There were no repetitive sequences detected in /x/temp_ins.fa\n"
        self.assertEqual(self.run_in_tmp(text, ['INS_A'])['INS_A'], [])
        text = HEADER + " 2000 10.0 0.0 0.0 INS_Z 1 300 (0) + AluSx SINE/Alu 2 301 (11) 1\n"
        self.assertEqual(self.run_in_tmp(text, ['INS_A'])['INS_A'], [])


if __name__ == '__main__':
    unittest.main()
