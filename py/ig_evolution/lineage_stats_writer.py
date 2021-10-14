import sys

from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

import utils
import dataset
import vj_annotator

class ClonalLineageStatWriter:
    def __init__(self, clonal_lineages):
        self.clonal_lineages = clonal_lineages
        self.dataset = self.clonal_lineages[0].Dataset()

    def _GetHighestMultiplicity(self, lineage):
        return max([self.dataset.GetSeqMultiplicity(seq.id) for seq in lineage.FullLengthSeqIdIter()])

    def _GetNumberNonTrivialSequences(self, lineage):
        return len([seq for seq in lineage.FullLengthSeqIdIter() if self.dataset.GetSeqMultiplicity(seq.id) > 1])

    def _GetRootSeq(self, lineage, abundant_v):
        num_shms_all = sys.maxint
        root_seq_all = SeqRecord(Seq(''), id = '')
        num_shms_prod = sys.maxint
        root_seq_prod = SeqRecord(Seq(''), id = '')
        for seq in lineage.FullLengthSeqIdIter():
            if utils.GetBaseGeneName(self.dataset.GetGeneHitBySeqName(seq.id, dataset.AnnotatedGene.V)) != abundant_v:
                continue
            cur_num_shms = len(self.dataset.GetVSHMsOutsideCDR3(seq.id)) + len(self.dataset.GetJSHMsOutsideCDR3(seq.id))
            if cur_num_shms < num_shms_all:
                num_shms_all = cur_num_shms
                root_seq_all = seq
            cur_seq = seq.seq #self.dataset.GetCDR3BySeqName(seq.id)
            aa_seq = str(Seq(cur_seq).translate())
            if aa_seq.find('*') != -1:
                continue
            if cur_num_shms < num_shms_prod:
                num_shms_prod = cur_num_shms
                root_seq_prod = seq
        if root_seq_prod.id != '':
            return root_seq_prod
        return root_seq_all

    def OutputStats(self, output_fname):
        fh = open(output_fname, 'w')
        fh.write('LineageID\tLineageSizeBeforeCleaning\tNumNonTrivialSeqs\tMaxMultiplicity\tClosestV\tClosestJ\tRootId\tRootSeq\tRootCDR3\tRootDistanceFromGermline\n')
        for l in sorted(self.clonal_lineages, key = lambda s : len(s), reverse = True):
            vj_ann = vj_annotator.VJGeneAnnotator(l)
            abundant_v = utils.GetBaseGeneName(vj_ann.GetAbundantGene(dataset.AnnotatedGene.V))
            abundant_j = utils.GetBaseGeneName(vj_ann.GetAbundantGene(dataset.AnnotatedGene.J))
            root_seq = self._GetRootSeq(l, abundant_v)
            num_shms_in_root = len(self.dataset.GetVSHMsOutsideCDR3(root_seq.id)) + len(self.dataset.GetJSHMsOutsideCDR3(root_seq.id))
            fh.write(l.id() + '\t' + str(len(l)) + '\t' + str(self._GetNumberNonTrivialSequences(l)) + '\t' + str(self._GetHighestMultiplicity(l)) + '\t' + abundant_v + '\t' + abundant_j + '\t' + root_seq.id + '\t' + root_seq.seq + '\t' + self.dataset.GetCDR3BySeqName(root_seq.id) + '\t' + str(num_shms_in_root) + '\n')
        fh.close()

    def OutputClonalLineageAssignment(self, output_fname):
        fh = open(output_fname, 'w')
        fh.write('SequenceID\tLineageID\n')
        for lineage in sorted(self.clonal_lineages, key = lambda s : len(s), reverse = True):
            for seq in lineage.FullLengthSeqIdIter():
                fh.write(str(lineage.id()) + '\t' + seq.id + '\n')
        fh.close()
