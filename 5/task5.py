"""
Piotr Dobrowolski
(291528)
Task 5
"""

from Bio import Seq
from Bio import SeqIO
from Bio.Blast import NCBIWWW,NCBIXML

import urllib
import urllib2

from BeautifulSoup import BeautifulStoneSoup 

import random
import time
import sys


minimal_blast_score = 0.0

class FastaLoader(object):
    """ Simply load sequences from fasta file """
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
        """ returns random sequence """
        return self.sequences.values()\
            [random.randint(0, len(self.sequences)-1)]

class PFAMManager(object):
    """ Manage every connections with PFAM """
    def __init__(self):
        self.url = "http://pfam.sanger.ac.uk/search/sequence"

    def parse_job(self, xml):
        """ From xml returned from PFAM parse job,
        I need only url, and only url is returning """
        soup = BeautifulStoneSoup(xml)
        try:
            urlname = soup.result_url.string
        except Exception as nt:
            return None
        return urlname

    def load_pfam(self, sequence):
        """ from sequence produce xml with results, returned from
        PFAM web page """
        #REQUEST A JOB
        post = {}
        post['seq'] = sequence
        #post['evalue'] = '10.0'
        #post['ga'] = '1'
        #post['searchBs'] = '0'
        post['output'] = 'xml'
        data = urllib.urlencode(post)
        request_job = urllib2.Request(self.url, data)
        try:
            response_job = urllib2.urlopen(request_job)
            xml_job_info = response_job.read()
        except urllib2.HTTPError as e:
            print e.reason
            return None
        except urllib2.URLError as e:
            print e.reason
            return None

        result_url = self.parse_job(xml_job_info)
        
        if not result_url:
            return None

        #WAIT FOR DONE JOB
        while(True):
            request_result = urllib2.Request(result_url)
            try:
                response_result = urllib2.urlopen(request_result)
                xml_result = response_result.read()
            except urllib2.HTTPError as e:
                print e.reason
                break
            except urllib2.URLError as e:
                print e.reason
                break

            if not (xml_result == ''):
                break
            else:
                """ Job is undone on the PFAM server, wait a second
                and try again. 
                WARNING: If HTTPError, URLError or returns XML then
                infinite loop occurred """
                time.sleep(1)

        return xml_result

    def parse_match(self, match):
        """ Parse single match from PFAM """
        if not match:
            return None
        return (match.attrs[0][1], match.attrs, match.seq.contents[1])

    def parse_result(self, xml):
        """ Parse all matches from PFAM """
        try:
            soup = BeautifulStoneSoup(xml)
        except Exception as e:
            return []
        matches = soup.findAll('match')
        return map(self.parse_match, matches)

    def search(self, sequence):
        """ Search all matches in PFAM,
        first load xml from server, than parse it"""
        xml = self.load_pfam(sequence)
        return self.parse_result(xml)

class ProteinMaker(object):
    """ Makes several protein sequences from single one """
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
    """ Manages connection to blast server using Bio module Blast """
    def __init__(self, sequence):
        self.seq = sequence

    def blast2subjects(self, blastrecs, minscore):
        """ Interpreting blastrec as list of subject, only with score greater
        than minscore """
        sbjcts = [] #pair (subject, score)
        for blastrec in blastrecs:
            for align in blastrec.alignments:
                for hsp in align.hsps:
                    """ add High Scoring Segment Pairs subject and score """
                    sbjcts.append((hsp.sbjct, hsp.score))
        return [sbjct for sbjct, score in sbjcts if  score >= minscore]

    def search(self):
        res = NCBIWWW.qblast("blastx", "nr", self.seq)
        blast_records = NCBIXML.parse(res)
        ali = self.blast2subjects(blast_records, minimal_blast_score)
        seqs = []
        for a in ali:
            newseq = Seq.Seq(a, Seq.Alphabet.ProteinAlphabet)
            seqs.append(newseq)
            
        return seqs


def find_function(fasta_file):
    """ Main Function, should does thinks like in task description """
    fl = FastaLoader(fasta_file)
    seq = fl.random_sequence()
    print "Chosen sequence: "
    print seq

    pfam_m = PFAMManager()
    found_using_blast = {}
    found_using_frames = {}
    
    print "Makes protein..."
    pm = ProteinMaker(seq)
    protein_seqs_frames = pm.make()
    print "Loading from BLAST..."
    bm = BLASTManager(seq)
    protein_seqs_blast = bm.search()

    print "Searching pfam, part1: frames..."
    lenn = len(protein_seqs_frames)
    i = 0
    for ps in protein_seqs_frames:
        if ps:
            found_all = pfam_m.search(ps)
            for ident, attr, seq in found_all:
                found_using_frames[ident] = (attr, seq)
        i += 1
        print "{0}%".format(100 * float(i)/float(lenn))

    print ''
    print "Searching pfam, part2: blast..."
    lenn = len(protein_seqs_blast)
    i = 0
    for ps in protein_seqs_blast:
        if ps:
            found_all = pfam_m.search(ps)
            for ident, attr, seq in found_all:
                found_using_blast[ident] = (attr, seq)
        i += 1
        print "{0}%".format(100 * float(i)/float(lenn))
    print ''
    
    all_sum = {}
    all_sum.update(found_using_frames)
    all_sum.update(found_using_blast)
    return (found_using_frames, found_using_blast,
        all_sum)

if __name__ == '__main__':
    """ Sequences from fasta files can be loaded via command line """
    for fn in sys.argv[1:]:
        print find_function(fn)        
