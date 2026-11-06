import tempfile,unittest
from pathlib import Path
from chimera_bench.tools.centrifuger_convert import convert_centrifuger_quant

class NativeProfileTests(unittest.TestCase):
    def convert(self,rows,method='native-sequence'):
        with tempfile.TemporaryDirectory() as tmp:
            src=Path(tmp)/'quant.tsv';out=Path(tmp)/'profile.tsv'
            raw='name\ttaxID\ttaxRank\tgenomeSize\tnumReads\tnumUniqueReads\tabundance\n'+rows
            src.write_text(raw)
            convert_centrifuger_quant(quant_path=src,out_path=out,profile_method=method)
            self.assertEqual(src.read_text(),raw)
            return {r.split('\t')[0]:float(r.split('\t')[-1]) for r in out.read_text().splitlines()[1:]}
    def test_native_length_and_positive_saved_abundance(self):
        got=self.convert('root\t1\troot\t100\t999\t0\t1\na\t2\tspecies\t100\t1\t0\t0.6\nb\t3\tspecies\t400\t999\t0\t0.4\nc\t4\tspecies\t100\t10\t0\t0\nstrain\t5\tstrain\t400\t999\t0\t0.4\n')
        self.assertEqual(set(got),{'2','3'});self.assertAlmostEqual(got['2'],100*60/220);self.assertAlmostEqual(sum(got.values()),100)
    def test_count_control_retains_zero_em_species(self):
        got=self.convert('a\t2\tspecies\t100\t1\t0\t1\nb\t3\tspecies\t400\t3\t0\t0\n','read-counts')
        self.assertEqual(got,{'2':25.,'3':75.})
    def test_reject_invalid_native_values(self):
        for length,abundance in [('0','0.1'),('-1','0.1'),('nan','0.1'),('100','nan'),('100','-0.1'),('100','inf')]:
            with self.subTest(length=length,abundance=abundance),self.assertRaises(ValueError):self.convert(f'a\t2\tspecies\t{length}\t1\t0\t{abundance}\n')
    def test_duplicate_rejected(self):
        with self.assertRaises(ValueError):self.convert('a\t2\tspecies\t100\t1\t0\t1\na\t2\tspecies\t100\t1\t0\t1\n')
    def test_empty_profile(self):
        self.assertEqual(self.convert(''),{})
    def test_malformed_rejected(self):
        with self.assertRaises(ValueError):self.convert('bad\trow\n')

if __name__=='__main__':unittest.main()
