'''
Created on Jul 26, 2017

@author: Marta Luksza, mluksza@ias.edu
'''
from Bio.Blast import NCBIXML
from Bio import pairwise2
from Bio.SubsMat import MatrixInfo as matlist
from math import log, exp

class Aligner:
    """
    Class to align neoantigens with IEDB epitopes and compute TCR-recognition
    probabilities.
    """
    INF = float("inf")
        
    @staticmethod
    def align(seq1, seq2):
        """
        Smith-Waterman alignment with default parameters.
        """
        matrix = matlist.blosum62
        gap_open = -11
        gap_extend = -1
        aln = pairwise2.align.localds(seq1.upper(), seq2.upper(), matrix, gap_open, gap_extend)
        return aln
    
    @staticmethod
    def log_sum(v):
        """
        compute the logarithm of a sum of exponentials
        """
        if len(v) == 0:
            return -Aligner.INF
        ma = max(v)
        if ma == -Aligner.INF:
            return -Aligner.INF
        return log(sum(map(lambda x: exp(x-ma), v))) + ma


    def __init__(self):
        # dictionary of computed Ri-values mapped to neoantigen identifiers
        self.Ri = {}
        # dictionary of IEDB epitope alignments mapped to neoantigen identifiers
        self.alignments = {}
        # dictionary of the highest scoring alignments mapped to neoantigen identifiers
        self.maximum_alignment = {}

    def read_all_blast_alignments(self, xm_lpath):
        """
        Read precomputed blastp alignments from xml files,
        compute alignment scores,
        find the highest scoring alignment for each neoantigen.
        """
        maxscore_by_epitope_id = dict()
        with open(xm_lpath) as f:
            blast_records = NCBIXML.parse(f)
            for brecord in blast_records:
                epitope_id = int(str(brecord.query).split("_")[1])
                for alignment in brecord.alignments:
                    if epitope_id not in self.alignments:
                        self.alignments[epitope_id] = {}
                        self.maximum_alignment[epitope_id] = 0
                        maxscore_by_epitope_id[epitope_id] = 0

                    species = " ".join((str(alignment).split())[1:-3])
                    for hsp in alignment.hsps:
                        if "-" not in hsp.query and "-" not in hsp.sbjct:
                            al = Aligner.align(hsp.query, hsp.sbjct)
                            if len(al) > 0:
                                al = al[0]
                                self.alignments[epitope_id][species] = al
                                if epitope_id not in maxscore_by_epitope_id:
                                    breakpoint()
                                # assert epitope_id in maxscore, maxscore
                                if al[2] > maxscore_by_epitope_id[epitope_id]:
                                    self.maximum_alignment[epitope_id] = species
                                    maxscore_by_epitope_id[epitope_id] = al[2]

    def compute_rval(self, a=26, k=1):
        """
        Compute TCR-recognition probabilities for each neoantigen.
        """
        # iterate over all neoantigens
        for i in self.alignments: 
            # energies of all bound states of neoantigen i
            binding_energies = list(map(lambda el: -k*(a-el[2]), self.alignments[i].values()))
            # partition function, over all bound states and an unbound state
            lZ = Aligner.log_sum(binding_energies + [0])
            lGb = Aligner.log_sum(binding_energies)
            R = exp(lGb-lZ)
            self.Ri[i] = R

    def get_rval(self, i):
        """
        Return precomputed R value and the highest scoring alignment
        for a given neoantigen i.
        """
        empty_alignment = [None, None, 0]
        if i in self.Ri:
            species = self.maximum_alignment[i]
            al = self.alignments[i][species]
            species = str(species).replace(" ", "_")
            return [self.Ri[i], species, al]
        return [0., None, empty_alignment]


    