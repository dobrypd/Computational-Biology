#!/usr/bin/env python

#Piotr Dobrowolski (291528)

#task3

from Bio import AlignIO

import urllib2

class TreeWrapper(object):
    def __init__(self):
        pass

    def load_from_pfam_seq(self, multiple_seq_alignment):
        pass

    def get_tree(self):
        pass

    def __str__(self):
        pass

class ProteinFamily(object):
    tree = TreeWrapper()

    def __load_from_pfam(self, identifier):
        pfam_http_addr_prefix = \
            'http://pfam.sanger.ac.uk/family/alignment/download/format?acc='
        pfam_http_addr_sufix = \
            '&alnType=seed&format=stockholm&order=t&case=l&gaps=dashes'
        self.identifier = identifier
        url = pfam_http_addr_prefix + identifier + pfam_http_addr_sufix
        http_handler = urllib2.urlopen(url)
        format = 'stockholm'
        self.multiple_seq_alignment = AlignIO.read(
                http_handler, format)

    def __init__(self, identifier):
        self.__load_from_pfam(identifier)
        self.tree.load_from_pfam_seq(self.multiple_seq_alignment)

    def as_MultipleSeqAlignment(self):
        return self.multiple_seq_alignment

    def get_tree(self):
        return self.tree.get_tree()

class TreeChanger(object):
    protein_family = None
    sub_matrix = None
    penalty = None
    
    def __init__(self, protein_faminy, substitution_matrix, penalty):
        pass

    def __next_cut_iteration(self):
        pass

    def improve(self, iterations):
        for i in iterations:
            print "Iteration {0}".format(i)
            self.__next_cut_iteration()
            print str(protein_family.get_tree())
