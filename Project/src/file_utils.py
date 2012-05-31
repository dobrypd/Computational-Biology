'''
Computational Biology
Main Project

file_utils.py - classes and methods operating on files.

Author: Piotr Dobrowolski
no:     291528
'''
from Bio import SeqIO


class FastaFileWrapper(object):
    ''' imply load sequences from fasta file.'''
    def __init__(self, filename):
        self.sequences = {}
        self._parse_file(self._load_file(filename))

    def _load_file(self, filename):
        '''Opens file, and catch exceptions.'''
        try:
            return open(filename)
        except IOError as(errno, strerr):
            print "I/O error while loading file({0}): {1}"\
            .format(errno, strerr)

    def _parse_file(self, file_fasta):
        '''Loading sequences into dictionary from file, closing file.'''
        fasta_sequences = SeqIO.parse(file_fasta, "fasta")
        for sequence in fasta_sequences:
            self.sequences[sequence.id] = sequence.seq
        file_fasta.close()

    def seuqences_gen(self):
        for _, v in self.sequences.items():
            yield v


def sequences_gen(filename):
    '''Returns generator with sequences from some (filename) fasta file.'''
    ffw = FastaFileWrapper(filename)
    return ffw.seuqences_gen()
