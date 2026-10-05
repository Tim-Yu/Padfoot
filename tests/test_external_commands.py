"""The external commands receive each argument as one word, whatever spaces or shell characters it holds (run: python3 -m
unittest discover -s tests). Stub RepeatMasker / minimap2 / samtools / bedtools executables record the argv bash hands
them."""
import os, sys, tempfile, unittest
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from padfoot import annot

STUB = r"""#!/bin/bash
line=$(basename "$0"); for a in "$@"; do line+=$'\t'"$a"; done
printf '%s\n' "$line" >> "$PADFOOT_ARGV_LOG"   # one write: minimap2 and samtools sort log concurrently
case "$(basename "$0") $1" in
  "samtools faidx") printf '>%s\nACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACG\n' "$3" ;;
  "samtools sort") cat > /dev/null ;;
  "bedtools intersect") cat "$3" > "$PADFOOT_ARGV_LOG.a" ;;   # what bedtools found in the -a file
esac
exit 0
"""


def make_sv(vcf_id, sv_type, bp_1, bp_2, ins_seq=''):
    return annot.SV(bp_1, '+', bp_2, '-', 10, '', sv_type, 0.4, vcf_id, False, False, ins_seq, '', '', 0, 0)


class ExternalCommandArguments(unittest.TestCase):
    def setUp(self):
        self.cwd = os.getcwd()
        self.d = tempfile.TemporaryDirectory()
        self.addCleanup(self.d.cleanup)
        self.addCleanup(os.chdir, self.cwd)
        bindir = os.path.join(self.d.name, 'bin')
        os.makedirs(bindir)
        for tool in ('RepeatMasker', 'minimap2', 'samtools', 'bedtools'):
            path = os.path.join(bindir, tool)
            with open(path, 'w') as fh:
                fh.write(STUB)
            os.chmod(path, 0o755)
        self.log = os.path.join(self.d.name, 'argv.log')
        env = {'PATH': bindir + os.pathsep + os.environ.get('PATH', ''), 'PADFOOT_ARGV_LOG': self.log}
        for key, value in env.items():
            self.addCleanup(self.restore_env, key, os.environ.get(key))
            os.environ[key] = value
        self.work = os.path.join(self.d.name, 'work dir')
        os.makedirs(self.work)
        os.chdir(self.work)

    @staticmethod
    def restore_env(key, value):
        if value is None:
            os.environ.pop(key, None)
        else:
            os.environ[key] = value

    def calls(self, tool):
        with open(self.log) as fh:
            argvs = [l.rstrip('\n').split('\t') for l in fh]
        return [a[1:] for a in argvs if a[0] == tool]

    def test_species_with_a_space_is_one_repeatmasker_argument(self):
        sv = make_sv('INS1', 'INS', ('chr19', 100), ('chr19', 100), ins_seq='ACGT' * 20)
        ref = os.path.join(self.d.name, 'ref dir', 'chr19.fa')
        for specie in ('Mus musculus', 'mus_musculus', 'human; touch pwned'):
            if os.path.exists(self.log):
                os.remove(self.log)
            self.assertTrue(annot.write_ins([sv], ref, 2, specie, True))
            self.assertEqual(self.calls('RepeatMasker'), [['-species', specie, 'temp_ins.fa']])
            self.assertEqual(self.calls('minimap2'),
                             [['-ax', 'map-ont', ref, 'temp_ins.fa', '-k', '17', '-y', '-K', '5G', '-t', '2', '--eqx']])
            self.assertEqual(self.calls('samtools'), [['sort', '-@', '2', '-m', '4G'], ['index', '-@', '2', 'temp_ins.bam']])
        self.assertFalse(os.path.exists('pwned'))   # the ';' was not run as a shell command

    def test_samtools_faidx_reference_and_region_are_single_arguments(self):
        ref = os.path.join(self.d.name, 'ref dir', 'chr19.fa')
        annot.get_microhomology([make_sv('DEL1', 'DEL', ('chr19', 1215001), ('chr19', 1223000))], ref)
        self.assertEqual(self.calls('samtools'), [['faidx', ref, 'chr19:1214976-1215026'],
                                                  ['faidx', ref, 'chr19:1222975-1223025']])

    def test_bedtools_gets_the_repeat_file_as_one_argument(self):
        sv = make_sv('DEL1', 'DEL', ('chr19', 1215001), ('chr19', 1223000))
        rm_bed = os.path.join(self.d.name, 'my annotations', 'rm.bed')
        annot.annot_bp_repeat({'DEL1': sv}, rm_bed)
        self.assertEqual(self.calls('bedtools'), [['intersect', '-a', 'temp_bps.bed', '-b', rm_bed, '-wb']])
        with open(self.log + '.a') as fh:   # temp_bps.bed was complete when bedtools ran (it used to be empty)
            self.assertEqual(fh.read(), 'chr19\t1214996\t1215006\tDEL1\tBP1\nchr19\t1222995\t1223005\tDEL1\tBP2\n')


if __name__ == '__main__':
    unittest.main()
