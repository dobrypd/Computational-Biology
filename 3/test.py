import unittest
import Bio
from Bio import AlignIO
from Bio import Phylo


from task3 import ProteinFamily
from task3 import TreeWrapper

#in test directory should be
identifier = 'PF05371'
format_ = 'stockholm'


class TestProteinFamily(unittest.TestCase):
    def setUp(self):
        pass
    def test_loading_from_pfam(self):

        pf = ProteinFamily(identifier)

        should_be_align = AlignIO.read(identifier + "." + format_, format_)
        should_be_align.sort()
        
        what_have = pf.as_MultipleSeqAlignment()
        what_have.sort()

        for i in xrange(len(should_be_align)):
            for j in xrange(len(should_be_align[i])):
                self.assertSequenceEqual(should_be_align[i][j],
                    what_have[i][j])

class TestTreeWrapper(unittest.TestCase):
    def setUp(self):
        self.seq = AlignIO.read(identifier + "." + format_, format_)

    def test_load_from_pfam_seq(self):
        tree = TreeWrapper()
        tree.load_from_pfam_seq(self.seq)
        phylo_tree = tree.get_tree()
        self.assertIsNotNone(phylo_tree, "No tree returned")
        self.assertTrue(phylo_tree.rooted, "Silly is tree test")


if __name__ == '__main__':
    unittest.main()
