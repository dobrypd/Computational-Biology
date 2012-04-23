"""
Piotr Dobrowolski
(291528)
Task 4
"""

import MarkovModel
#from Bio.HMM import MarkovModel
from Bio.Alphabet import Alphabet
from Bio import SeqIO

from itertools import product
from math import floor

""" size of emissions """
size = 4

class StatesAlphabet(Alphabet):
    promoter = 'p'
    out_promoter = 'o'
    letters = [
        promoter,
        out_promoter
    ]

class EmissionAlphabet(Alphabet):
    letters = []
    def __init__(self, k):
        nucl = ['A', 'C', 'T', 'G']
        self.letters = ["".join(i) for i in product(nucl, repeat=k)]


def get_avg_promoter_len(promoters):
    sum_ = sum([len(i) for i in promoters])
    avg = sum_ / len(promoters)
    return avg

def iterate_other_alphabet(seq, k):
    """ changes to k-th len letters alphabet """
    for i in xrange(int(floor(len(seq)/k))):
        yield seq[i*k:(i+1)*k]

def count_emissions_from_promoters(emissions, promoters, letter_size):
    """ count no of emission for all possible emissions """
    ret = {}
    for emission in emissions:
        ret[emission] = 0

    for promoter in promoters:
        for e in iterate_other_alphabet(promoter, letter_size):
            ret[e] += 1

    return ret

def count_emissions_out_promoters(emissions, genome, letter_size):
    """ count no of emission for all possible emissions """
    ret = {}
    for emission in emissions:
        ret[emission] = 0

    for e in iterate_other_alphabet(genome, letter_size):
        ret[e] += 1

    return ret

def create_HMM(promoters_file, genome_file):
    """ 
    Creates hidden markov model
    promoters_file - fasta file handler with promoters
    genome_file - fasta file handler with sequences of all genome
    returns HMM
    """
    letters_in_alp_len = size
    states = StatesAlphabet()
    emissions = EmissionAlphabet(letters_in_alp_len)

    hmm_builder = MarkovModel.MarkovModelBuilder(states, emissions)
    hmm_builder.allow_all_transitions()
    hmm_builder.set_random_probabilities()
    
    promoters_seqs = SeqIO.parse(promoters_file, "fasta")
    genome_seq = SeqIO.read(genome_file, "fasta")
    
    promoters = [i.seq.tostring() for i in promoters_seqs]
    genome = genome_seq.seq.tostring()

    #SET TRANSITINS SCORES
    ######################
    #promoter -> no promoter
    hmm_builder.set_transition_score(states.promoter, states.out_promoter,
            float(letters_in_alp_len) / float(get_avg_promoter_len(promoters)))
    #no promoter -> promoter 
    hmm_builder.set_transition_score(states.out_promoter, states.promoter,
            float(len(promoters)) * letters_in_alp_len / len(genome))
    #stay in this same state
    hmm_builder.set_transition_score(states.promoter, states.promoter,
            1.0-float(letters_in_alp_len) / float(
                get_avg_promoter_len(promoters)))
    hmm_builder.set_transition_score(states.out_promoter, states.out_promoter,
            1.0-float(letters_in_alp_len) / float(
                get_avg_promoter_len(promoters)))
    
    #SET EMISSIONS SCORES
    #####################
    #EMISSIONS FROM PROMOTER
    promoters_emission_dict = count_emissions_from_promoters(
            emissions.letters, promoters, letters_in_alp_len)
    emission_all_in_promoters = sum(len(i) for i in promoters)
    for emission, count in promoters_emission_dict.items():
        hmm_builder.set_emission_score(states.promoter, 
                emission, float(count) / emission_all_in_promoters)
    #EMISSIONS OUTSIDE PROMOTER
    no_promoters_emission_dict = count_emissions_out_promoters(
            emissions.letters, genome, letters_in_alp_len)
    emissions_all_out_promoters = len(genome) / letters_in_alp_len
    for emission, count in no_promoters_emission_dict.items():
        hmm_builder.set_emission_score(states.out_promoter,
                emission, float(count) / emissions_all_out_promoters)

    return hmm_builder.get_markov_model()
    

def promoters_position(hmm, genome_file):
    """
    Using hidden markov model, predict positions of promoters
    returns list of positions
    """
    genome_seq = SeqIO.read(genome_file, "fasta")
    genome = genome_seq.seq.tostring()

    states_seq = hmm.viterbi(
            [i for i in iterate_other_alphabet(genome, size)],
            StatesAlphabet())
    all_pos =  [i*size for i, x in enumerate(states_seq[0]) \
        if x == StatesAlphabet().promoter]

    new_pos = []
    last = -1
    for i in all_pos:
        if not i == last:
            last = i
            new_pos.append(i)
        last += 1
    return new_pos

if __name__ == '__main__':

    import pdb
    import sys
    if not len(sys.argv[1:]) == 3:
        print """ open with args :
            promoters_sequences_fasta_file 
            complite_genome_fasta_file 
            find_in_fasta_file
        """
        sys.exit(1)

    file_prom = sys.argv[1]
    file_gen_0 = sys.argv[2]
    file_gen_1 = sys.argv[3]

    print "creating hmm..."
    hmm = create_HMM(file_prom, file_gen_0)
    pdb.set_trace()
    print "promoters position searching..."
    pos_list = promoters_position(hmm, file_gen_1)
    print pos_list
    pdb.set_trace()

