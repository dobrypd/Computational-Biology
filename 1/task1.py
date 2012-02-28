#!/usr/bin/env python


from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import SingleLetterAlphabet


class ProbesGenerator(object):
    '''
    Class Probes Generator 
    (generates probes for current file in fasta)
    '''
    # file with sequences in fasta format
    filename = ''
    # reverse complement DNA sequences, and probes dictionaries (ID is key)
    sequences = {}
    probes = {}


    def __init__(self, filename):
        '''
        ProbesGenerator constructor
        trying to open and parse input file
        '''
        self.filename = filename
        self.__parse_file(self.__load_file())

    def __load_file(self):
        '''
        Opening file, and exception catch
        '''
        try:
            return open(self.filename)
        except IOError as (errno, strerr):
            print "I/O error while loading file({0}): {1}"\
                .format(errno, strerr)

    def __parse_file(self, file_fasta):
        '''
        Loading sequences into dictionary from file, closing file
        '''
        fasta_sequences = SeqIO.parse(file_fasta, "fasta")
        for sequence in fasta_sequences:
            self.sequences[sequence.id] = sequence.seq.\
                                          reverse_complement().tostring()
        file_fasta.close()

    def __create_subsets(self, lenght):
        subsets = {}
        for name, sequence in self.sequences.items():
            subsets[name] = set()
            for i in (sequence[i:i+lenght] \
                    for i in xrange(len(sequence) - lenght + 1)):
                subsets[name].add(i)
        return subsets

    def __try_to_create_probes(self, subsets, f):
        for name, sequence in self.sequences.items():
            sets_difference = subsets[name].copy()
            for name2, sequence2 in self.sequences.items():
                if name != name2:
                    sets_difference -= subsets[name2]
            #check if difference is not empty (probe is correct)
            if len(sets_difference) > 0:
                # chceck function f if not None
                if (f != None):
                    could_be_added = set()
                    for probe in sets_difference:
                        if f(probe):
                            could_be_added.add(probe)
                    if len(could_be_added) > 0:
                        self.probes[name] = could_be_added
                else:
                     self.probes[name] = sets_difference

    def __raise_cannot_be_found(self):
        all_sequences_ids = set(self.sequences.keys())
        correct_probes_ids = set(self.probes.keys())
        cannot_be_distinguish = all_sequences_ids - correct_probes_ids;
        raise Exception("Cannot distinguish differences in sequences: %s" \
                    % [i for i in cannot_be_distinguish])

    def __prepare_seq_records_list(self):
            ret = []
            for probe_id, probe_seq_as_str_list in self.probes.items():
                ret.append(SeqRecord(Seq(probe_seq_as_str_list.pop(), \
                                SingleLetterAlphabet()), \
                            id=probe_id, name=probe_id, \
                            description=("Probe for sequence %s" % (probe_id))))
            return ret
    
    def get_probes(self, k=None, f=None):
        '''
        get_probes
            k - lenght of the probes, if None then try to find smallest
            f - function, check if seq is a good probe
            
        '''
        self.probes = {}
        found = False
        lenght = 1
        if k:
            lenght = k

        #if cannot be found with min_len then pass 
        #(because every probes has to be this same lenght)
        min_len = min( (len(x) for x in self.sequences.values()) )

        while (not found) and (lenght <= min_len):
            #subset is dictionary (key is ID) 
            #   (value is list of sequence subset lenght `lenght`)
            subsets = {}
            subsets = self.__create_subsets(lenght)
            self.__try_to_create_probes(subsets, f)

            # check if all sequences has it's own probe
            found = True
            for name, sequence in self.sequences.items():
                if not self.probes.has_key(name):
                    found = False
            lenght = lenght + 1
            if k:
                break
        
       
        # check if all sequences has it's own probe
        if not found:
            self.__raise_cannot_be_found()
        else:
            return self.__prepare_seq_records_list()



def get_probes(filename, k=None, f=None):
    '''
    get_probes generates
        filename - in FASTA format
        k - lenght of the probes, if None then try to find smallest
        f - function, check if seq is a good probe
        
    '''
    PG = ProbesGenerator(filename)
    return PG.get_probes(k, f)
