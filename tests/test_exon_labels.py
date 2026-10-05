"""Unit tests for the exon/intron/flank labels, the in-frame test and the known-fusion lookup
(run: python3 -m unittest discover -s tests)."""
import os, sys, unittest
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from padfoot import annot

# a 12-exon gene on chr1:1,000-23,000, exons of 1 kb every 2 kb (exon numbers in transcript order)
def gene(strand):
    starts = [1000 + 2000 * i for i in range(12)]
    ends = [s + 999 for s in starts]
    nums = list(range(1, 13)) if strand == '+' else list(range(12, 0, -1))
    return [starts, ends, ['exon%d' % k for k in nums], [strand]]

class ExonLabels(unittest.TestCase):
    def test_plus_strand(self):
        g = gene('+')
        self.assertEqual(annot.exon_label(g, 500, '+'), 'promoter/utr')     # 5' flank
        self.assertEqual(annot.exon_label(g, 1500, '+'), 'exon1')
        self.assertEqual(annot.exon_label(g, 2500, '+'), 'intron1')
        self.assertEqual(annot.exon_label(g, 20500, '+'), 'intron10')       # not 'intron0'
        self.assertEqual(annot.exon_label(g, 22500, '+'), 'intron11')
        self.assertEqual(annot.exon_label(g, 23500, '+'), 'exon12')
        self.assertEqual(annot.exon_label(g, 25000, '+'), 'downstream')     # 3' flank, was 'intron2'
    def test_minus_strand_numbers_follow_the_transcript(self):
        g = gene('-')
        self.assertEqual(annot.exon_label(g, 500, '-'), 'downstream')       # coordinate-first = 3' end
        self.assertEqual(annot.exon_label(g, 2500, '-'), 'intron11')        # between exon12 and exon11
        self.assertEqual(annot.exon_label(g, 22500, '-'), 'intron1')        # between exon2 and exon1
        self.assertEqual(annot.exon_label(g, 25000, '-'), 'promoter/utr')   # coordinate-last = 5' end
    def test_single_exon_gene(self):
        g = [[1000], [2000], ['exon1'], ['+']]
        self.assertEqual(annot.exon_label(g, 1500, '+'), 'exon1')
        self.assertEqual(annot.exon_label(g, 900, '+'), 'promoter/utr')
        self.assertEqual(annot.exon_label(g, 2100, '+'), 'downstream')

class FakeSV:
    def __init__(self, genes, typ):
        self.genes = [genes[0], genes[1], typ]; self.score = 3; self.impact = 'HIGH'; self.cancer = [' ', ' ']

class KnownFusions(unittest.TestCase):
    def test_partner_listed_under_either_gene(self):
        fusion = {'BCR': ['ABL1', 'JAK2']}
        for pair in (('BCR', 'ABL1'), ('ABL1', 'BCR')):
            sv = FakeSV([(pair[0], 'intron1', '+'), (pair[1], 'intron1', '+')], 'possible_fusion')
            annot.cancer_annot_svs([sv], {}, fusion)
            self.assertEqual((sv.genes[2], sv.impact, sv.score), ('oncogenic_fusion', 'HIGH_ONCO', 6), pair)
    def test_unlisted_pair_stays_a_possible_fusion(self):
        sv = FakeSV([('BCR', 'intron1', '+'), ('TP53', 'intron1', '+')], 'possible_fusion')
        annot.cancer_annot_svs([sv], {}, {'BCR': ['ABL1']})
        self.assertEqual((sv.genes[2], sv.score), ('possible_fusion', 3))

if __name__ == '__main__':
    unittest.main()
