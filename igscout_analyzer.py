import os
import sys
import shutil

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.pairwise2 import format_alignment
from Bio import pairwise2

from enum import Enum

import matplotlib as mplt
mplt.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

import seaborn as sns
import numpy as np

###################################################
def PrepareOutputDir(output_dir):
    if os.path.exists(output_dir):
        shutil.rmtree(output_dir)
    os.mkdir(output_dir)

def ReadFasta(fasta_fname):
    records = []
    for r in SeqIO.parse(fasta_fname, 'fasta'):
        r.seq = str(r.seq).upper()
        records.append(r)
    print str(len(records)) + ' records were extracted from ' + fasta_fname
    return records

def GetBaseGeneName(d_name):
    return d_name.split('*')[0]

###################################################
class SegmentType(Enum):
        UNKNOWN_TYPE = 0
        KNOWN_GENE = 1
        AMBIGUOUS_ALLELE = 2
        AMBIGUOUS_GENE = 3
        NOVEL_VARIATION = 4
        NOVEL_GENE = 5
        INACCURATE_GENE = 6

class DSegmentClassification:
    def __init__(self, segment_type, segment_alignment, gene_alignment, gene_ids):
        self.segment_type = segment_type
        self.segment_alignment = segment_alignment
        self.gene_alignment = gene_alignment
        self.gene_ids = gene_ids
        self.segment = ''.join([s for s in segment_alignment if s != '-'])

    def GetDs(self):
        if isinstance(self.gene_ids, (list,)):
            return ','.join([d_id for d_id in self.gene_ids])
        return self.gene_ids

    def GetSegment(self):
        return self.segment
    
    def GetSegmentType(self):
        return self.segment_type

    def __str__(self):
        classification_str = str(self.segment_type) + ', closest D genes: ' + self.GetDs() + '\n'
        classification_str += self.segment_alignment + '\n'
        classification_str += self.gene_alignment
        return classification_str

class DSegmentClassifier:
    def __init__(self, d_genes):
        self.d_genes = d_genes
        self.relative_diff_pos = 0.2
        self.max_variation_diff = 0.25

    def _GetContainingDGenes(self, d_segment):
        d_indices = []
        for i in range(len(self.d_genes)):
            if self.d_genes[i].seq.find(d_segment) != -1:
                d_indices.append(i)
        return d_indices

    def _AlignNovelSegment(self, d_segment):
        best_alignment = ''
        best_score = 0
        best_id = ''
        for d in self.d_genes:
            alignment = pairwise2.align.localms(d_segment, d.seq, 5, -4, -3, -.1)[0]
            if alignment[2] > best_score:
                best_score = alignment[2]
                best_alignment = alignment
                best_id = d.id
        return best_alignment, best_id

    def _ComputeNonGappedRange(self, segment_gene_alignment):
        start_pos = -1
        for i in range(len(segment_gene_alignment[0])):
            start_pos = i
            if segment_gene_alignment[0][i] != '-' and segment_gene_alignment[1][i] != '-':
                break
        end_pos = -1
        for i in range(len(segment_gene_alignment[0])):
            end_pos = len(segment_gene_alignment[0]) - i - 1
            if segment_gene_alignment[0][end_pos] != '-' and segment_gene_alignment[1][end_pos] != '-':
                break
        return (start_pos, end_pos + 1)

    def _ClassifyNovelSegment(self, segment_gene_alignment):
        non_gapped_range = self._ComputeNonGappedRange(segment_gene_alignment)
        inner_diff = []
        for i in range(non_gapped_range[0], non_gapped_range[1]):
            if segment_gene_alignment[0][i] != segment_gene_alignment[1][i]:
                inner_diff.append(i)
        if len(inner_diff) == 0:
            return SegmentType.INACCURATE_GENE
        #for pos in inner_diff:
        # TODO: add checking how close diff positions to start / end
        if float(len(inner_diff)) / (non_gapped_range[1] - non_gapped_range[0]) > self.max_variation_diff:
            return SegmentType.NOVEL_GENE
        return SegmentType.NOVEL_VARIATION

    def _SegmentIsAmbiguousGene(self, d_indices):
        d_base_names = set([GetBaseGeneName(self.d_genes[ind].id) for ind in d_indices])
        return len(d_base_names) > 1

    def _ClassifySubstringSegment(self, segment, d_indices):
        if self._SegmentIsAmbiguousGene(d_indices):
            return SegmentType.AMBIGUOUS_GENE
        if len(d_indices) > 1:
            return SegmentType.AMBIGUOUS_ALLELE
        return SegmentType.KNOWN_GENE

    def _GetAlignmentSeqForSubstring(self, substring, string):
        substr_index = string.index(substring)
        return '-' * substr_index + substring + '-' * (len(string) - len(substring) - substr_index)

    def _CreateSubstringSegmentClassification(self, segment_type, segment, d_indices):
        if segment_type == SegmentType.AMBIGUOUS_GENE:
            return DSegmentClassification(segment_type, segment, '', [GetBaseGeneName(self.d_genes[ind].id) for ind in d_indices])
        main_d_gene = self.d_genes[d_indices[0]]
        segment_alignment_seq = self._GetAlignmentSeqForSubstring(segment, main_d_gene.seq)
        return DSegmentClassification(segment_type, segment_alignment_seq, main_d_gene.seq, main_d_gene.id) 

    def Classify(self, d_segment):
        d_indices = self._GetContainingDGenes(d_segment)
#        print d_indices
        if len(d_indices) == 0:
            d_hit_alignment, d_hit = self._AlignNovelSegment(d_segment)
            return DSegmentClassification(self._ClassifyNovelSegment(d_hit_alignment), d_hit_alignment[0], d_hit_alignment[1], d_hit)
        segment_type = self._ClassifySubstringSegment(d_segment, d_indices)
        return self._CreateSubstringSegmentClassification(segment_type, d_segment, d_indices)

def OutputClassificationsToDF(classifications, output_fname):
    fh = open(output_fname, 'w')
    fh.write('Segment\tType\tClosest_Ds\n')
    for c in classifications:
        fh.write(c.GetSegment() + '\t' + str(c.GetSegmentType()) + '\t' + c.GetDs() + '\n')
    fh.close()

def main(segment_fasta, d_gene_fasta, output_dir):
    PrepareOutputDir(output_dir)
    segments = ReadFasta(segment_fasta)
    d_genes = ReadFasta(d_gene_fasta)
    segment_classifier = DSegmentClassifier(d_genes)
    classifications = []
    for s in segments:
        classification = segment_classifier.Classify(s.seq)
        print str(classification) + '\n'
        classifications.append(classification)
    OutputClassificationsToDF(classifications, os.path.join(output_dir, 'segment_annotation.txt'))

if __name__ == '__main__':
    if len(sys.argv) != 4:
        print "Invalid input parameters"
        print "python igscout_analyzer.py inferred_segments.fasta IGHD.fasta output_dir"
        sys.exit(1)
    main(sys.argv[1], sys.argv[2], sys.argv[3])
