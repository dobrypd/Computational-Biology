"""
Piotr Dobrowolski
(291528)
Task 5
"""

from Bio import SeqIO
from Bio.Blast import NCBIWWW,NCBIXML

import random
import sys


minimal_blast_score = 350.0

class FastaLoader(object):
    def __init__(self, filename):
        random.seed()
        self.sequences = {}
        self._parse_file(self._load_file(filename))

    def _load_file(self, filename):
        '''
        Opening file, and exception catch
        '''
        try:
            return open(filename)
        except IOError as (errno, strerr):
            print "I/O error while loading file({0}): {1}"\
                .format(errno, strerr)

    def _parse_file(self, file_fasta):
        '''
        Loading sequences into dictionary from file, closing file
        '''
        fasta_sequences = SeqIO.parse(file_fasta, "fasta")
        for sequence in fasta_sequences:
            self.sequences[sequence.id] = sequence.seq
        file_fasta.close()

    def random_sequence(self):
        return self.sequences.values()\
            [random.randint(0, len(self.sequences)-1)]

class PFAMManager(object):
    def __init__(self):
        pass

    def search(self, sequence):
        return []

class ProteinMaker(object):
    def __init__(self, sequence):
        self.seq = sequence
    def make(self):
        """Prepare all possible proteins"""
        normal_seq = self.seq
        reversed_seq = self.seq.complement()
        return [normal_seq.translate(),
                normal_seq[1:].translate(),
                normal_seq[2:].translate(),
                reversed_seq.translate(),
                reversed_seq[1:].translate(),
                reversed_seq[2:].translate()]
                
class BLASTManager(object):
    def __init__(self, sequence):
        self.seq = sequence
    def blast2subjects(self, blastrec, minscore):
        """ Interpreting blastrec as list of subject, only with score greater
        than minscore """
        sbjcts = [] #pair (subject, score)
        for align in blastrec.alignments:
            for hsp in align.hsps:
                """ add High Scoring Segment Pairs subject and score """
                sbjcts.append((hsp.sbjct, hsp.score))
        return [sbjct for sbjct, score in sbjcts if  score >= minscore]

    def search(self):
        res = NCBIWWW.qblast("blastn","nr",self.seq)
        blast_records = NCBIXML.parse(res)
        return self.blast2subjects(blast_records, minimal_blast_score)


def find_function(fasta_file):
    fl = FastaLoader(fasta_file)
    seq = fl.random_sequence()
    print "Chosen sequence: "
    print seq

    pfam_m = PFAMManager()
    found_using_blast = set()
    found_using_frames = set()
    
    pm = ProteinMaker(seq)
    protein_seqs_frames = pm.make()
    bm = BLASTManager(seq)
    protein_seqs_blast = bm.search()

    for ps in protein_seqs_frames:
        found_using_frames.add(
                pfam_m.search(ps))

    for ps in protein_seqs_blast:
        found_using_blast.add(
                pfam_m.search(ps))
        
    return (found_using_frames, found_using_frames,
            found_using_blast | found_using_frames)

if __name__ == '__main__':
    for fn in sys.argv[1:]:
        print find_function(fn)        
