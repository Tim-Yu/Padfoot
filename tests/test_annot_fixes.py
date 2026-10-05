"""Unit tests for pad_unlisted_segments() and get_repeat() (run: python3 -m unittest discover -s tests)."""
import logging, os, shutil, sys, tempfile, unittest
from collections import defaultdict
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from padfoot import annot

logger = logging.getLogger()


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

    def get_cna(self, records, svs):
        with tempfile.TemporaryDirectory() as d:
            vcf = os.path.join(d, 'w.vcf')
            with open(vcf, 'w') as fh:
                fh.write(WAKHAN_HEADER)
                fh.writelines(records)
            with self.assertLogs(level='WARNING') as logs:
                cnas, base = annot.get_CNA(vcf, svs, ploidy=2.0)
                logger.warning('sentinel')   # assertLogs needs at least one record
        return cnas, base, [r.getMessage() for r in logs.records if r.getMessage() != 'sentinel']

    def test_header_only_vcf_is_copy_neutral(self):
        # Wakhan writes altered segments only, so a copy-neutral tumour gives a header-only VCF; this used to crash in
        # int(np.median([])) (ValueError: cannot convert float NaN to integer)
        sv = make_sv(10)
        cnas, base, warnings = self.get_cna([], [sv])
        self.assertEqual(base, [1, 1])
        for hp in (1, 2):
            self.assertEqual([(c.pos_1, c.pos_2, c.cn, c.dir1) for c in cnas[('chr3', hp)]],
                             [(1, 198295559, 1.0, 'NEUT')])
        self.assertFalse(sv.cn_altering)
        self.assertEqual(len(warnings), 1, warnings)
        self.assertIn('lists no copy-number segment', warnings[0])

    def test_per_copy_coverage_below_one_does_not_divide_by_zero(self):
        # COV1/CN1 = 0.5 -> int(median) * 0.75 = 0 -> round(sv.supp / 0) used to raise ZeroDivisionError
        sv = make_sv(10)
        record = wakhan_record('chr3', 1000001, 2000000, 1.0, 1.0).replace(':15.0:15.0\n', ':0.5:0.5\n')
        cnas, _, warnings = self.get_cna([record], [sv])
        self.assertEqual([(c.pos_1, c.pos_2, c.cn) for c in cnas[('chr3', 1)]],
                         [(1, 1000000, 1.0), (1000001, 2000000, 1.0), (2000001, 198295559, 1.0)])
        self.assertFalse(sv.cn_altering)
        self.assertEqual(len(warnings), 1, warnings)
        self.assertIn('per-copy coverage', warnings[0])

    def test_cn_altering_estimate_is_unchanged_with_normal_coverage(self):
        sv = make_sv(22)
        _, _, warnings = self.get_cna([wakhan_record('chr3', 1000001, 2000000, 2.0, 1.0)], [sv])
        self.assertEqual(sv.cn_altering, 2)    # round(22 / (int(30 / 2) * 0.75))
        self.assertEqual(warnings, [])


def make_sv(supp):
    return annot.SV(('chr3', 5000), '+', ('chr3', 9000), '-', supp, '', 'DEL', 0.4, 'DEL1', False, False, '', '', '', 0, 0)


class ContigLengthsFromFai(unittest.TestCase):
    """A Wakhan VCF without ##contig lines gives pysam contigs without a length; the padding then takes the lengths
    from the reference .fai and must give exactly what the ##contig lines give."""
    CONTIGS = '##contig=<ID=chr3,length=198295559>\n'
    FAI = 'chr3\t198295559\t6\t60\t61\nchr4\t190214555\t201600330\t60\t61\nchrM\t16569\t394997123\t60\t61\n'
    RECORDS = [wakhan_record('chr3', 87000001, 94000000, 0.0, 0.0), wakhan_record('chr3', 162807955, 162908546, 2.0, 1.0)]
    GENES = {'chr3': [['BEFORE', 'IN_LOSS', 'IN_GAIN', 'AFTER'], [1000000, 90000000, 162850000, 190000000],
                      [1050000, 90100000, 162860000, 190050000]],
             'chr4': [['NO_SEGMENT'], [5000000], [5100000]]}

    def run_cna(self, with_contigs, records, fai):
        header = WAKHAN_HEADER if with_contigs else WAKHAN_HEADER.replace(self.CONTIGS, '')
        if with_contigs:
            header = header.replace(self.CONTIGS, self.CONTIGS + '##contig=<ID=chr4,length=190214555>\n')
        with tempfile.TemporaryDirectory() as d:
            vcf = os.path.join(d, 'w.vcf')
            with open(vcf, 'w') as fh:
                fh.write(header)
                fh.writelines(records)
            ref_fai = None
            if fai:
                ref_fai = os.path.join(d, 'ref.fa.fai')
                with open(ref_fai, 'w') as fh:
                    fh.write(self.FAI)
            with self.assertLogs(level='INFO') as logs:
                cnas, base = annot.get_CNA(vcf, [], ploidy=2.0, ref_fai=ref_fai)
        by_gene = {}
        annot.annot_CNAs(self.GENES, cnas, base, by_gene)
        segs = {k: [(c.pos_1, c.pos_2, c.cn) for c in v] for k, v in cnas.items()}
        genes = {g: (tuple(b.CN), tuple(b.CN_impact)) for g, b in by_gene.items()}
        warnings = [r.getMessage() for r in logs.records if r.levelname == 'WARNING']
        return segs, genes, warnings

    def test_missing_contig_lines_are_padded_from_the_fai_exactly_as_with_them(self):
        segs_ref, genes_ref, warn_ref = self.run_cna(True, self.RECORDS, fai=False)
        segs, genes, warnings = self.run_cna(False, self.RECORDS, fai=True)
        self.assertEqual(segs, segs_ref)
        self.assertEqual(genes, genes_ref)
        self.assertEqual(genes['BEFORE'], ((1.0, 1.0), ('NEUT', 'NEUT')))     # not the chr3 loss wrapped round (DEL/DEL)
        self.assertEqual(genes['IN_LOSS'], ((0.0, 0.0), ('DEL', 'DEL')))
        self.assertEqual(genes['AFTER'], ((1.0, 1.0), ('NEUT', 'NEUT')))      # not NA beyond the last segment
        self.assertEqual(genes['NO_SEGMENT'], ((1.0, 1.0), ('NEUT', 'NEUT')))  # chromosome without a record
        self.assertNotIn(('chrM', 1), segs)                                    # only primary or listed contigs are padded
        self.assertEqual(warn_ref, [])
        self.assertEqual(len(warnings), 1, warnings)
        self.assertIn('ref.fa.fai', warnings[0])

    def test_header_only_vcf_without_contig_lines(self):
        segs, genes, warnings = self.run_cna(False, [], fai=True)
        self.assertEqual(segs[('chr3', 1)], [(1, 198295559, 1.0)])
        self.assertEqual(segs[('chr4', 2)], [(1, 190214555, 1.0)])
        self.assertEqual(set(genes.values()), {((1.0, 1.0), ('NEUT', 'NEUT'))})

    def test_no_length_anywhere_is_reported(self):
        _, genes, warnings = self.run_cna(False, self.RECORDS, fai=False)
        self.assertTrue(any('No length for contig(s) chr3' in w for w in warnings), warnings)


class AnnotCNAs(unittest.TestCase):
    def cnas_for(self, records, contigs):
        header = WAKHAN_HEADER.replace('##contig=<ID=chr3,length=198295559>\n',
                                       ''.join(f'##contig=<ID={c},length={l}>\n' for c, l in contigs))
        with tempfile.TemporaryDirectory() as d:
            vcf = os.path.join(d, 'w.vcf')
            with open(vcf, 'w') as fh:
                fh.write(header)
                for r in records:
                    fh.write(wakhan_record(*r))
            return annot.get_CNA(vcf, [])

    def test_par_gene_keeps_the_chrx_copy_number(self):
        cnas, ploidy = self.cnas_for([('chrX', 1000001, 2000000, 2.0, 1.0), ('chrY', 1, 57227415, 0.0, 0.0)],
                                     [('chrX', 156040895), ('chrY', 57227415)])
        genes = {'chrX': [['SHOX'], [614344], [669411]], 'chrY': [['SHOX'], [614344], [669411]]}
        by_gene = {}
        annot.annot_CNAs(genes, cnas, ploidy, by_gene)
        self.assertEqual(by_gene['SHOX'].ref_id, 'chrX')
        self.assertEqual(by_gene['SHOX'].CN, [1.0, 1.0])          # neutral on chrX, not the chrY 0/0

    def test_gene_beyond_the_contig_end_is_left_unset(self):
        cnas, ploidy = self.cnas_for([('chr19', 5000001, 6000000, 2.0, 1.0)], [('chr19', 58617616)])
        genes = {'chr19': [['PARD6G', 'STK11'], [80147232, 1167558], [80257514, 1238431]]}   # PARD6G is really on chr18
        by_gene = {}
        annot.annot_CNAs(genes, cnas, ploidy, by_gene)
        self.assertEqual(by_gene['PARD6G'].CN, ['NA', 'NA'])       # untouched default: copy number unknown (printed as NA)
        self.assertEqual(by_gene['STK11'].CN, [1.0, 1.0])          # neutral, outside the 5-6 Mb gain


@unittest.skipUnless(shutil.which('bedtools'), 'bedtools not in PATH')
class AnnotBpRepeat(unittest.TestCase):
    def test_breakpoints_in_repeats_are_annotated(self):
        # temp_bps.bed was still open (unflushed) when bedtools read it, so no breakpoint ever got a repeat class
        svs = {'DEL1': make_sv(10), 'INS1': annot.SV(('chr3', 20000), '+', ('chr3', 20000), '-', 5, '', 'INS', 0.4, 'INS1',
                                                       False, False, 'ACGT', '', '', 0, 0)}
        cwd = os.getcwd()
        with tempfile.TemporaryDirectory() as d:
            os.chdir(d)
            try:
                with open('rm.bed', 'w') as fh:
                    fh.write('chr3\t4000\t4500\tLINE/L1\n'        # no breakpoint
                             'chr3\t8990\t9300\tSINE/Alu\n'       # DEL1 BP2 (9000)
                             'chr3\t19000\t20003\tLTR/ERVL\n')    # INS1 (20000)
                annot.annot_bp_repeat(svs, os.path.join(d, 'rm.bed'))
            finally:
                os.chdir(cwd)
        self.assertEqual(svs['DEL1'].repeat_bp, ['', 'SINE/Alu'])
        self.assertEqual(svs['INS1'].repeat_bp, ['LTR/ERVL', ''])


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
