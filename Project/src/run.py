'''
Computational Biology
Main Project

run.py - loading files, executing other modules

Author: Piotr Dobrowolski
no:     291528
'''

import sys


def print_help():
    print """Usage: {0} NEW_ORGRANISM GENOME POSITIONS FUNCTIONS
    where
NEW_ORGANISK: fasta file with DNA of organism to analyze
GENOME: genome sequence, fasta file
POSITIONS: genome positions, .bed file
FUNCTIONS: genome annotations (genome->functions), .go file"""\
    .format(sys.argv[0])

if __name__ == "__main__":
    if not len(sys.argv) == 5:
        print_help()
        exit(1)
    for i in sys.argv[1:]:
        print i
