import unittest

from task5 import find_function, FastaLoader, BLASTManager, ProteinMaker

testing_seq_fn="./testing_seq.fa"

class FindFunctionTest(unittest.TestCase):
    def x_test_find_run_with_file(self):
        self.assertIsNotNone(find_function(testing_seq_fn))

class TestBlastManager(unittest.TestCase):
    def test_slow_not_mocked(self):
        fl = FastaLoader(testing_seq_fn)
        bm = BLASTManager(fl.random_sequence())
        print bm.search()

    def test_pm(self):
        fl = FastaLoader(testing_seq_fn)
        pm = ProteinMaker(fl.random_sequence())
        print pm.make()


if __name__ == '__main__':
    unittest.main()
