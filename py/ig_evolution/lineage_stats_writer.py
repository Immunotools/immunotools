import sys
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
        num_shms = sys.maxint
        root_seq = ''
        for seq in lineage.FullLengthSeqIdIter():
            if utils.GetBaseGeneName(self.dataset.GetGeneHitBySeqName(seq.id, dataset.AnnotatedGene.V)) != abundant_v:
                continue
            cur_num_shms = len(self.dataset.GetVSHMsOutsideCDR3(seq.id)) + len(self.dataset.GetJSHMsOutsideCDR3(seq.id))
            if cur_num_shms < num_shms:
                num_shms = cur_num_shms
                root_seq = seq
        return root_seq

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
