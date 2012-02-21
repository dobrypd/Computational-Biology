#! /usr/bin/env python

from Bio import SeqIO

"""
get_probes
    fileName - nazwa pliku w formacie FASTA
    k - dlugosc sond, None - najmniejsza
    f - funkcja, ktora sprawdza wiezy zwracanych sond
        funkcja dla sekwencji sprawdza czy moze wystapic w rozwiazaniu
"""
def get_probes(fileName, k=None, f=None):
    file = open(fileName)

    seqsF = SeqIO.parse(file, "fasta")
    seqs = {}
    fseq = {} #to return
    for seq in seqsF:
        seqs[seq.id] = seq.seq.reverse_complement() #TODO: check

    #znalezc unikalne podslowa
    #i zrobic odwrotnie komplementarne .reverse_complement()
    found = 0
    lenght = 1
    while found == 0:
        cutSets = {};
        for name, seq in seqs.items():
            #tworze zbiory wszystkich podslow
            #i odejmuje
            cutSets[name] = set()
            #wszystkie podslowa dl k
            #[s[i:i+k] for i in len(s)-k] for k in 

            for i in range (len(seq) - lenght + 1):
                cutSets[name].add(seq[i:i+(lenght-1)])

        #generate:
        for name, seq in seqs.items():
            if not fseq.has_key(name):
                ret = cutSets[name]
                for name2, seq2 in seqs.items():
                    if name != name2:
                        ret -= cutSets[name2]
                    if len(ret) > 0:
                        fseq[name] = ret
        found = 1
        for name, s in seqs.items():
            if not fseqs.has_key(name):
                found = 0

        lenght = lenght + 1
#    if ...
#    else:a
#        raise Exception("Cannot be done")
 
    for name, r in fseq.items():
        print name
        print r

    

