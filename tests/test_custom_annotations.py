"""Unit tests for the custom --gff / --rm annotation inputs (run: python3 -m unittest discover -s tests)."""
import gzip, os, sys, tempfile, unittest
from unittest import mock
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from padfoot import annot, preprocess
from padfoot import main as padfoot_main

GENCODE_LIKE = """##gff-version 3
chr1\tHAVANA\tgene\t100\t900\t.\t+\t.\tID=g1;gene_type=protein_coding;gene_name=GENE_A
chr1\tHAVANA\ttranscript\t100\t900\t.\t+\t.\tID=t1;Parent=g1;gene_type=protein_coding;gene_name=GENE_A;transcript_name=GENE_A-201
chr1\tHAVANA\texon\t100\t300\t.\t+\t.\tID=e1;Parent=t1;gene_type=protein_coding;gene_name=GENE_A;transcript_name=GENE_A-201;exon_number=1
chr1\tHAVANA\tCDS\t150\t300\t.\t+\t0\tID=c1;Parent=t1;gene_type=protein_coding;gene_name=GENE_A;transcript_name=GENE_A-201;exon_number=1
chr1\tHAVANA\texon\t600\t900\t.\t+\t.\tID=e2;Parent=t1;gene_type=protein_coding;gene_name=GENE_A;transcript_name=GENE_A-201;exon_number=2
chr1\tHAVANA\tgene\t2000\t3000\t.\t+\t.\tID=g9;gene_type=lncRNA;gene_name=LNC_X
chr1\tHAVANA\ttranscript\t2000\t3000\t.\t+\t.\tID=t9;Parent=g9;gene_type=lncRNA;gene_name=LNC_X;transcript_name=LNC_X-201
chr1\tHAVANA\texon\t2000\t3000\t.\t+\t.\tID=e9;Parent=t9;gene_type=lncRNA;gene_name=LNC_X;transcript_name=LNC_X-201;exon_number=1
chr2\tHAVANA\tgene\t10\t50\t.\t-\t.\tID=g3;gene_type=protein_coding;gene_name=GENE_C
chr2\tHAVANA\ttranscript\t10\t50\t.\t-\t.\tID=t3;Parent=g3;gene_type=protein_coding;gene_name=GENE_C;transcript_name=GENE_C-201
chr2\tHAVANA\texon\t10\t50\t.\t-\t.\tID=e4;Parent=t3;gene_type=protein_coding;gene_name=GENE_C;transcript_name=GENE_C-201;exon_number=1
"""
EXPECTED_TABLE = [['chr1', 'GENE_A', '+', 'gene', '100', '900'], ['chr1', 'GENE_A', '+', 'GENE_A-201', '100', '900'],
                  ['chr1', 'GENE_A', '+', 'exon1', '100', '300'], ['chr1', 'GENE_A', '+', 'exon2', '600', '900'],
                  ['chr2', 'GENE_C', '-', 'gene', '10', '50'], ['chr2', 'GENE_C', '-', 'GENE_C-201', '10', '50'],
                  ['chr2', 'GENE_C', '-', 'exon1', '10', '50']]

# a real RepeatMasker .out: 3 header lines, right-aligned columns, a complement hit and the optional overlap flag '*'
RM_OUT = """   SW   perc perc perc  query      position in query           matching       repeat              position in  repeat
score   div. del. ins.  sequence    begin     end    (left)    repeat         class/family         begin  end (left)   ID

  463    1.3  0.6  1.7  chr1          10001     10468 (248945954) +  (TAACCC)n      Simple_repeat            1    471    (0)      1
 3612   11.4 21.5  1.3  chr1          10469     11447 (248944975) C  TAR1           Satellite/telo       (399)   1712    483      2
  484   25.1 13.2  0.0  chr1          11505     11675 (248944747) C  L1MC5a         LINE/L1             (2382)   5648   5452      3 *
"""
RM_BED = [['chr1', '10001', '10468', 'Simple_repeat'], ['chr1', '10469', '11447', 'Satellite/telo'],
          ['chr1', '11505', '11675', 'LINE/L1']]


def write(path, text, gz=False):
    with (gzip.open if gz else open)(path, 'wt') as fh:
        fh.write(text)
    return path


def rows(path):
    with preprocess.open_text(path) as fh:
        return [l.rstrip('\n').split('\t') for l in fh]


class GenerateGff(unittest.TestCase):
    def setUp(self):
        self.d = tempfile.TemporaryDirectory()
        self.addCleanup(self.d.cleanup)

    def path(self, name):
        return os.path.join(self.d.name, name)

    def test_gzipped_gff3_is_read_and_the_gz_output_is_what_get_genes_reads(self):
        src = write(self.path('in.gff3.gz'), GENCODE_LIKE, gz=True)
        out = self.path('gff_file.gff3.gz')
        self.assertEqual(preprocess.generate_gff(src, out), 2)        # the lncRNA is dropped
        with open(out, 'rb') as fh:
            self.assertEqual(fh.read(2), b'\x1f\x8b')                # gzip, like beds/hg38.gff3.gz
        self.assertEqual(rows(out), EXPECTED_TABLE)
        genes, exon_pos = annot.get_genes(out)
        self.assertEqual(genes['chr1'][0], ['GENE_A'])
        self.assertEqual(genes['chr2'][0], ['GENE_C'])
        self.assertEqual(exon_pos['GENE_A'][2], ['exon1', 'exon2'])

    def test_plain_and_gzipped_inputs_give_the_same_table(self):
        out_gz, out_plain = self.path('a.gff3.gz'), self.path('b.gff3.gz')
        preprocess.generate_gff(write(self.path('in.gff3.gz'), GENCODE_LIKE, gz=True), out_gz)
        preprocess.generate_gff(write(self.path('in.gff3'), GENCODE_LIKE), out_plain)
        self.assertEqual(rows(out_gz), rows(out_plain))

    def test_get_genes_also_reads_a_plain_table(self):
        out = self.path('table.gff3')                                  # no .gz suffix -> plain text
        preprocess.generate_gff(write(self.path('in.gff3'), GENCODE_LIKE), out)
        with open(out, 'rb') as fh:
            self.assertNotEqual(fh.read(2), b'\x1f\x8b')
        self.assertEqual(annot.get_genes(out)[0]['chr1'][0], ['GENE_A'])

    def test_non_gencode_gff_yields_no_gene(self):
        refseq = "chr1\tBestRefSeq\tgene\t100\t900\t.\t+\t.\tID=gene-A;gene=A;gene_biotype=protein_coding\n"
        self.assertEqual(preprocess.generate_gff(write(self.path('r.gff3'), refseq), self.path('r.out.gz')), 0)


class GenerateRm(unittest.TestCase):
    def setUp(self):
        self.d = tempfile.TemporaryDirectory()
        self.addCleanup(self.d.cleanup)

    def run_rm(self, text, gz=False):
        src = write(os.path.join(self.d.name, 'in.gz' if gz else 'in.txt'), text, gz=gz)
        out = os.path.join(self.d.name, 'rm.bed')
        n = preprocess.generate_rm(src, out)
        result = rows(out)
        self.assertEqual(n, len(result))
        return result

    def test_repeatmasker_out(self):
        self.assertEqual(self.run_rm(RM_OUT), RM_BED)

    def test_gzipped_repeatmasker_out(self):
        self.assertEqual(self.run_rm(RM_OUT, gz=True), RM_BED)

    def test_tab_separated_out_rows_without_header_still_work(self):
        # the only layout the old pandas reader accepted
        text = ''.join('\t'.join(l.split()) + '\n' for l in RM_OUT.splitlines()[3:])
        self.assertEqual(self.run_rm(text), RM_BED)

    def test_bed4_and_wider_bed(self):
        bed4 = ''.join('\t'.join(r) + '\n' for r in RM_BED)
        self.assertEqual(self.run_rm(bed4), RM_BED)
        bed6 = 'track name=rmsk\n' + ''.join('\t'.join(r + ['0', '+']) + '\n' for r in RM_BED)
        self.assertEqual(self.run_rm(bed6), RM_BED)    # extra columns dropped: annot_bp_repeat reports the LAST column

    def test_unrecognised_row_is_an_error(self):
        with self.assertRaises(ValueError):
            self.run_rm("chr1 start end Simple_repeat\n")


class MainResolvesRelativeAnnotationPaths(unittest.TestCase):
    """padfoot.main.main() chdirs into <out>/temp before reading --gff / --rm; relative paths must still resolve."""

    def setUp(self):
        self.cwd = os.getcwd()
        self.d = tempfile.TemporaryDirectory()
        self.addCleanup(self.d.cleanup)
        self.addCleanup(os.chdir, self.cwd)
        os.chdir(self.d.name)
        for name in ('sv.vcf', 'cna.vcf', 'ref.fa', 'ref.fa.fai'):
            write(name, '')
        write('custom.gff3.gz', GENCODE_LIKE, gz=True)
        write('custom.fa.out', RM_OUT)

    def run_main(self, *extra):
        captured = {}

        def annotate_things(args):
            captured['args'] = args
            return [], [], [], [], []
        argv = ['padfoot.py', '--sv-vcf', 'sv.vcf', '--cna-file', 'cna.vcf', '--ref', 'ref.fa', '--out-dir', 'out',
                '--genome', 'chm13', '--ploidy', '2', '--skip_RepeatMasker'] + list(extra)
        with mock.patch.object(sys, 'argv', argv), \
             mock.patch.object(padfoot_main, '_enable_logging'), \
             mock.patch.object(padfoot_main, '_check_external_dependencies'), \
             mock.patch.object(padfoot_main, 'annotate_things', side_effect=annotate_things):
            padfoot_main.main()
        return captured['args']

    def test_relative_gff_and_rm(self):
        args = self.run_main('--gff', 'custom.gff3.gz', '--rm', 'custom.fa.out')
        out = os.path.join(self.d.name, 'out')
        self.assertEqual(os.path.realpath(args.new_gff), os.path.realpath(os.path.join(self.d.name, 'custom.gff3.gz')))
        self.assertEqual(os.path.realpath(args.gff_file), os.path.realpath(os.path.join(out, 'gff_file.gff3.gz')))
        self.assertEqual(annot.get_genes(args.gff_file)[0]['chr1'][0], ['GENE_A'])
        self.assertEqual(os.path.realpath(args.rm_file), os.path.realpath(os.path.join(out, 'rm.bed')))
        self.assertEqual(rows(args.rm_file), RM_BED)

    def test_missing_gff_or_rm_is_reported(self):
        for extra, message in ((('--gff', 'missing.gff3'), 'GFF file does not exist'),
                               (('--rm', 'missing.out'), 'RepeatMasker file does not exist')):
            with self.assertRaises(SystemExit) as cm, self.assertLogs(level='ERROR') as logs:
                self.run_main(*extra)
            self.assertEqual(cm.exception.code, 1)
            self.assertIn(message, logs.output[0])
            os.chdir(self.d.name)

    def test_bundled_annotations_without_gff_or_rm(self):
        args = self.run_main('--genome', 'hg38')
        beds = os.path.realpath(os.path.join(os.path.dirname(__file__), '..', 'beds'))
        self.assertEqual(os.path.realpath(args.gff_file), os.path.join(beds, 'hg38.gff3.gz'))
        self.assertEqual(os.path.realpath(args.rm_file), os.path.join(beds, 'hg38_rm.bed.gz'))


if __name__ == '__main__':
    unittest.main()
