import os
import sys
import shutil
import operator
from Bio import SeqIO
import igraph 
from igraph import *

import dataset

class FullLengthClonalLineage:
    def __init__(self, cdr3_lineage):
        self.cdr3_lineage = cdr3_lineage
        self.dataset = cdr3_lineage.dataset

#    def _FindRoot(self):
#        self.root = -1
#        num_shms = sys.maxint
#        for v in self.VertexIterator():
#            seq = self.GetFullLengthSequenceByIndex(v)
#            cur_num_shms = 0
#            for gene_type in dataset.AnnotatedGene:
#                gene_shms = self.dataset.GetSHMsBySeqName(seq.id, gene_type)
#                cur_num_shms += len(gene_shms)
#            if cur_num_shms < num_shms:
#                num_shms = cur_num_shms
#                self.root = v

    def id(self):
        return self.cdr3_lineage.id()

    def CDR3Iter(self):
        return self.cdr3_lineage.CDR3Iter()

    def CDR3Length(self):
        for cdr3 in self.CDR3Iter():
            return len(cdr3)
        return -1

    def FullLengthSeqIdIter(self):
        return self.cdr3_lineage.FullLengthIter()

    def GetFullLengthSequenceByIndex(self, ind):
        return self.cdr3_lineage.GetFullLengthSequenceByIndex(ind)

    def GetFullLengthSequenceByName(self, seq_name):
        return self.dataset.GetSeqByName(seq_name)

    def GetCDR3ByFullLengthSequenceIndex(self, ind):
        fl_seq = self.GetFullLengthSequenceByIndex(ind)
        return self.dataset.GetCDR3BySeqName(fl_seq.id)

    def Dataset(self):
        return self.dataset

    def SHMDepth(self, seq_id):
        num_shms = 0
        for gene_type in dataset.AnnotatedGene:
            gene_shms = self.dataset.GetSHMsBySeqName(seq_id, gene_type)
            num_shms += len(gene_shms)
        return num_shms

    def __len__(self):
        return self.cdr3_lineage.NumFullLengthSequences()
