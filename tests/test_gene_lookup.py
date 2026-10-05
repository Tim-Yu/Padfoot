"""Unit tests for the breakpoint-to-gene lookup (run: python3 -m unittest discover -s tests)."""
import gzip, os, sys, unittest
from collections import defaultdict
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from padfoot import annot

BEDS = os.path.join(os.path.dirname(__file__), '..', 'beds')


class FakeSV:
    def __init__(self):
        self.genes = []


class SortedGeneLists(unittest.TestCase):
    def test_sort_gene_lists_orders_by_start_and_tracks_the_running_max_end(self):
        out = annot.sort_gene_lists([['B', 'A', 'C'], [500, 100, 300], [600, 900, 350]])
        self.assertEqual(out[0], ['A', 'C', 'B'])
        self.assertEqual(out[1], [100, 300, 500])
        self.assertEqual(out[2], [900, 350, 600])
        self.assertEqual(out[3], [900, 900, 900])

    def test_find_gene_at_resolves_nested_and_overlapping_genes(self):
        # padded spans (gene +/- 10 kb): BIG 100k-500k (body 110k-490k) contains SMALL 200k-250k (body 210k-240k);
        # LATE 480k-700k (body 490k-690k) overlaps BIG's end
        genels = annot.sort_gene_lists([['BIG', 'SMALL', 'LATE'], [100000, 200000, 480000], [500000, 250000, 700000]])
        names = genels[0]
        self.assertEqual(names[annot.find_gene_at(genels, 220000)], 'SMALL')   # innermost body
        self.assertEqual(names[annot.find_gene_at(genels, 245000)], 'BIG')     # SMALL's flank only, BIG's body wins
        self.assertEqual(names[annot.find_gene_at(genels, 300000)], 'BIG')     # intron of BIG after SMALL ended
        self.assertEqual(names[annot.find_gene_at(genels, 485000)], 'BIG')     # LATE's flank only, BIG's body wins
        self.assertEqual(names[annot.find_gene_at(genels, 495000)], 'LATE')    # LATE's body beats BIG's flank
        self.assertEqual(names[annot.find_gene_at(genels, 105000)], 'BIG')     # flank hit with no body candidate
        self.assertIsNone(annot.find_gene_at(genels, 99999))
        self.assertIsNone(annot.find_gene_at(genels, 700001))

    def test_bundled_gene_lists_are_sorted_for_every_chromosome(self):
        for name in ('hg38.gff3.gz', 'mm10.gff3.gz'):
            genes, _ = annot.get_genes(os.path.join(BEDS, name))
            for chrom, genels in genes.items():
                starts = genels[1]
                self.assertTrue(all(starts[i] <= starts[i + 1] for i in range(len(starts) - 1)), f"{name} {chrom}")
                self.assertEqual(len(genels), 4)


class RealBreakpoints(unittest.TestCase):
    """Cases from the P215003155 runs that the old two-bisect test got wrong or answered by luck."""

    @classmethod
    def setUpClass(cls):
        cls.genes, cls.exon_pos = annot.get_genes(os.path.join(BEDS, 'hg38.gff3.gz'))

    def _annotate(self, chrom, pos):
        sv = FakeSV()
        by_gene = defaultdict(list)
        annot.annotBPs(sv, (chrom, pos), self.genes, self.exon_pos, by_gene, 1)
        return sv.genes[0]

    def test_savana_id_74048_is_sppl2b_intron3_on_both_breakends(self):
        # was blank with the 25fd6d4 annotation: bisect_left on the unsorted end list landed on ENSG00000273734
        self.assertEqual(self._annotate('chr19', 2347430)[:2], ('SPPL2B', 'intron3'))
        self.assertEqual(self._annotate('chr19', 2347837)[:2], ('SPPL2B', 'intron3'))

    def test_severus_bnd6167_is_znf595_promoter(self):
        # upstream of chr4's first gene; blank while LMLN sat at the front of the chr4 list
        self.assertEqual(self._annotate('chr4', 89170)[:2], ('ZNF595', 'promoter/utr'))

    def test_relabelled_last_genes_are_found_on_their_own_chromosome(self):
        # the last gene of each chromosome used to be labelled with the next chromosome; probe each gene's midpoint
        expected = {'PARD6G': 'chr18', 'MZF1': 'chr19', 'RABL2B': 'chr22', 'FRG2B': 'chr10', 'PRMT2': 'chr21', 'PCMTD2': 'chr20'}
        rows = {}
        for line in gzip.open(os.path.join(BEDS, 'hg38.gff3.gz'), 'rt'):
            f = line.rstrip('\n').split('\t')
            if f[3] == 'gene' and f[1] in expected:
                rows[f[1]] = (f[0], (int(f[4]) + int(f[5])) // 2)
        for gene, chrom in expected.items():
            self.assertEqual(rows[gene][0], chrom, gene)
            self.assertEqual(self._annotate(chrom, rows[gene][1])[0], gene)

    def test_gene_body_beats_a_neighbours_flank(self):
        # synthetic fixture DEL chr19:1,215,001-1,223,000: the end lies in STK11 exon 10 and 5 kb before CBARP (within its
        # 10 kb flank); the exon wins. 1,100,500 lies in no gene body, so the innermost flank (SBNO2) is reported.
        self.assertEqual(self._annotate('chr19', 1223000)[:2], ('STK11', 'exon10'))
        self.assertEqual(self._annotate('chr19', 1215001)[:2], ('STK11', 'intron3'))
        self.assertEqual(self._annotate('chr19', 1100500)[:2], ('SBNO2', 'promoter/utr'))
        self.assertEqual(self._annotate('chr19', 900000)[:2], ('R3HDM4', 'intron6'))

    def test_position_outside_every_gene_is_blank(self):
        self.assertEqual(self._annotate('chr19', 20000), ())


if __name__ == '__main__':
    unittest.main()
