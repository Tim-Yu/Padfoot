"""Unit tests for the bundled data files and the annotation generator (run: python3 -m unittest discover -s tests)."""
import gzip, os, sys, tempfile, unittest
from collections import defaultdict
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from padfoot import annot, preprocess

BEDS = os.path.join(os.path.dirname(__file__), '..', 'beds')
HG38 = {'chr1': 248956422, 'chr2': 242193529, 'chr3': 198295559, 'chr4': 190214555, 'chr5': 181538259, 'chr6': 170805979,
        'chr7': 159345973, 'chr8': 145138636, 'chr9': 138394717, 'chr10': 133797422, 'chr11': 135086622, 'chr12': 133275309,
        'chr13': 114364328, 'chr14': 107043718, 'chr15': 101991189, 'chr16': 90338345, 'chr17': 83257441, 'chr18': 80373285,
        'chr19': 58617616, 'chr20': 64444167, 'chr21': 46709983, 'chr22': 50818468, 'chrX': 156040895, 'chrY': 57227415, 'chrM': 16569}
MM10 = {'chr1': 195471971, 'chr2': 182113224, 'chr3': 160039680, 'chr4': 156508116, 'chr5': 151834684, 'chr6': 149736546,
        'chr7': 145441459, 'chr8': 129401213, 'chr9': 124595110, 'chr10': 130694993, 'chr11': 122082543, 'chr12': 120129022,
        'chr13': 120421639, 'chr14': 124902244, 'chr15': 104043685, 'chr16': 98207768, 'chr17': 94987271, 'chr18': 90702639,
        'chr19': 61431566, 'chrX': 171031299, 'chrY': 91744698, 'chrM': 16299}


class CancerGeneTable(unittest.TestCase):
    def test_bundled_table_loads_without_network(self):
        self.assertTrue(os.path.exists(annot.CANCER_GENES_PATH))
        cancer_genes, fusion = annot.load_cancer_genes()
        self.assertIn('TSG', cancer_genes['TP53'])   # 'oncogene,TSG' in the table
        self.assertIn('oncogene', cancer_genes['MYC'])
        self.assertTrue(len(cancer_genes) > 500 and len(fusion) > 50)
        self.assertIn('KMT2A', fusion['ABI1'])

    def test_missing_columns_are_reported(self):
        fh = tempfile.NamedTemporaryFile('w', suffix='.tsv', delete=False)
        fh.write("Gene_symbol\tRole\nTP53\tTSG\n")
        fh.close()
        self.addCleanup(os.remove, fh.name)
        with self.assertRaises(ValueError):
            annot.load_cancer_genes(fh.name)


class BundledAnnotations(unittest.TestCase):
    def _check(self, name, lens):
        rows = [l.rstrip('\n').split('\t') for l in gzip.open(os.path.join(BEDS, name), 'rt')]
        beyond = [r for r in rows if r[0] in lens and int(r[5]) > lens[r[0]]]
        self.assertEqual(beyond, [], f"{name}: rows beyond their chromosome end, e.g. {beyond[:2]}")
        chroms_per_gene_block = defaultdict(set)
        block = None
        for r in rows:
            if r[3] == 'gene':
                block = (r[0], r[1], r[4])
            chroms_per_gene_block[block].add(r[0])
        self.assertTrue(all(len(c) == 1 for c in chroms_per_gene_block.values()), f"{name}: a gene block spans two chromosomes")
        return rows

    def test_hg38_rows_lie_within_their_chromosome(self):
        rows = self._check('hg38.gff3.gz', HG38)
        genes = {(r[0], r[1]) for r in rows if r[3] == 'gene'}
        self.assertIn(('chr18', 'PARD6G'), genes)   # used to be labelled chr19
        self.assertIn(('chrX', 'WASH6P'), genes)    # used to be labelled chrY
        self.assertIn(('chrY', 'WASH6P'), genes)    # used to be labelled chrM
        self.assertNotIn(('chr19', 'PARD6G'), genes)

    def test_mm10_rows_lie_within_their_chromosome(self):
        rows = self._check('mm10.gff3.gz', MM10)
        genes = {(r[0], r[1]) for r in rows if r[3] == 'gene'}
        self.assertIn(('chr1', 'Cr2'), genes)       # used to be labelled chr2
        self.assertIn(('chrY', 'Gm21748'), genes)   # used to be labelled chrM

    def test_get_genes_places_the_relabelled_genes_on_their_chromosome(self):
        genes, _ = annot.get_genes(os.path.join(BEDS, 'hg38.gff3.gz'))
        self.assertIn('PARD6G', genes['chr18'][0])
        self.assertNotIn('PARD6G', genes['chr19'][0])
        self.assertTrue(all(end <= HG38['chr19'] + 10000 for end in genes['chr19'][2]))


GENCODE_LIKE = """##gff-version 3
chr1\tHAVANA\tgene\t100\t900\t.\t+\t.\tID=g1;gene_type=protein_coding;gene_name=GENE_A
chr1\tHAVANA\ttranscript\t100\t900\t.\t+\t.\tID=t1;Parent=g1;gene_type=protein_coding;gene_name=GENE_A;transcript_name=GENE_A-201
chr1\tHAVANA\texon\t100\t300\t.\t+\t.\tID=e1;Parent=t1;gene_type=protein_coding;gene_name=GENE_A;transcript_name=GENE_A-201;exon_number=1
chr1\tHAVANA\texon\t600\t900\t.\t+\t.\tID=e2;Parent=t1;gene_type=protein_coding;gene_name=GENE_A;transcript_name=GENE_A-201;exon_number=2
chr1\tHAVANA\tgene\t5000\t6000\t.\t-\t.\tID=g2;gene_type=protein_coding;gene_name=GENE_B
chr1\tHAVANA\ttranscript\t5000\t6000\t.\t-\t.\tID=t2;Parent=g2;gene_type=protein_coding;gene_name=GENE_B;transcript_name=GENE_B-201
chr1\tHAVANA\texon\t5000\t6000\t.\t-\t.\tID=e3;Parent=t2;gene_type=protein_coding;gene_name=GENE_B;transcript_name=GENE_B-201;exon_number=1
chr2\tHAVANA\tgene\t10\t50\t.\t+\t.\tID=g3;gene_type=protein_coding;gene_name=GENE_C
chr2\tHAVANA\ttranscript\t10\t50\t.\t+\t.\tID=t3;Parent=g3;gene_type=protein_coding;gene_name=GENE_C;transcript_name=GENE_C-201
chr2\tHAVANA\texon\t10\t50\t.\t+\t.\tID=e4;Parent=t3;gene_type=protein_coding;gene_name=GENE_C;transcript_name=GENE_C-201;exon_number=1
"""


class GenerateGff(unittest.TestCase):
    def test_last_gene_of_a_chromosome_keeps_its_chromosome_and_the_last_gene_is_written(self):
        src = tempfile.NamedTemporaryFile('w', suffix='.gff3', delete=False)
        src.write(GENCODE_LIKE)
        src.close()
        out = src.name + '.padfoot.gff3'
        self.addCleanup(os.remove, src.name)
        self.addCleanup(lambda: os.path.exists(out) and os.remove(out))
        preprocess.generate_gff(src.name, out)
        rows = [l.rstrip('\n').split('\t') for l in open(out)]
        genes = {r[1]: r[0] for r in rows if r[3] == 'gene'}
        self.assertEqual(genes, {'GENE_A': 'chr1', 'GENE_B': 'chr1', 'GENE_C': 'chr2'})   # GENE_B was 'chr2', GENE_C was missing
        self.assertEqual([r[3] for r in rows if r[1] == 'GENE_A'], ['gene', 'GENE_A-201', 'exon1', 'exon2'])


if __name__ == '__main__':
    unittest.main()
