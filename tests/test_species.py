"""--specie handling: the mm10 hint and the log line when the cancer gene annotation is skipped (run: python3 -m unittest
discover -s tests)."""
import os, sys, tempfile, unittest
from types import SimpleNamespace
from unittest import mock
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from padfoot import annot
from padfoot import main as padfoot_main


class Mm10SpeciesHint(unittest.TestCase):
    def setUp(self):
        self.cwd = os.getcwd()
        self.d = tempfile.TemporaryDirectory()
        self.addCleanup(self.d.cleanup)
        self.addCleanup(os.chdir, self.cwd)
        os.chdir(self.d.name)
        for name in ('sv.vcf', 'cna.vcf', 'ref.fa', 'ref.fa.fai'):
            open(name, 'w').close()

    def run_main(self, specie):
        argv = ['padfoot.py', '--sv-vcf', 'sv.vcf', '--cna-file', 'cna.vcf', '--ref', 'ref.fa', '--out-dir', 'out',
                '--genome', 'mm10', '--specie', specie, '--ploidy', '2', '--skip_RepeatMasker']
        with mock.patch.object(sys, 'argv', argv), \
             mock.patch.object(padfoot_main, '_enable_logging'), \
             mock.patch.object(padfoot_main, '_check_external_dependencies'), \
             mock.patch.object(padfoot_main, 'annotate_things', return_value=([], [], [], [], [])):
            padfoot_main.main()
        os.chdir(self.d.name)

    def test_mouse_gets_the_repeatmasker_term_recommended(self):
        with self.assertLogs(level='WARNING') as logs:
            self.run_main('mouse')
        hint = [m for m in logs.output if 'mm10' in m]
        self.assertEqual(len(hint), 1, logs.output)
        self.assertIn('--specie mus_musculus', hint[0])

    def test_mus_musculus_spellings_need_no_hint(self):
        for specie in ('mus_musculus', 'Mus musculus'):
            with self.assertNoLogs(level='WARNING'):
                self.run_main(specie)


class CancerAnnotationSkipIsLogged(unittest.TestCase):
    def run_annotate(self, specie):
        args = SimpleNamespace(cna_vcf='cna.vcf', vcf_file='sv.vcf', threads=1, ref='ref.fa', out_dir='out', gff_file='g',
                               sv_caller='severus', cna_caller='wakhan', ploidy=2.0, specie=specie, rm_file='rm.bed',
                               run_repeatmasker=False, cancer_genes=None)
        stubs = {name: mock.DEFAULT for name in ('annot_CNAs', 'annot_SVS', 'annot_ins', 'get_microhomology',
                                                 'output_svs', 'output_genes', 'cancer_annot')}
        with mock.patch.multiple(annot, get_genes=mock.Mock(return_value=({}, {})), get_SVs=mock.Mock(return_value=[]),
                                 get_CNAs=mock.Mock(return_value=({}, [1, 1])), **stubs) as mocks:
            annot.annotate_things(args)
        return mocks['cancer_annot']

    def test_non_human_species_logs_the_skip(self):
        with self.assertLogs(level='WARNING') as logs:
            cancer_annot = self.run_annotate('mus_musculus')
        cancer_annot.assert_not_called()
        self.assertEqual(len(logs.output), 1, logs.output)
        self.assertIn('Cancer gene annotation skipped', logs.output[0])
        self.assertIn('mus_musculus', logs.output[0])

    def test_human_runs_the_cancer_annotation_without_a_warning(self):
        with self.assertNoLogs(level='WARNING'):
            cancer_annot = self.run_annotate('human')
        cancer_annot.assert_called_once()


if __name__ == '__main__':
    unittest.main()
