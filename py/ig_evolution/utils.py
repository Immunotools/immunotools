import os
import sys
import shutil
from enum import Enum

import matplotlib as mplt
mplt.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

from Bio import SeqIO
from Bio import pairwise2
from Bio.Seq import Seq

#########################################################################
def HammingDistance(seq1, seq2):
    if len(seq1) != len(seq2):
        print("ERROR: sequences " + seq1 + ' & ' + seq2 + ' have different lengths')
        sys.exit(1)
    dist = 0
    for i in range(len(seq1)):
        if seq1[i] != seq2[i]:
            dist += 1
    return dist

def EditDistance(seq1, seq2):
    alignment = pairwise2.align.globalxx(Seq(seq1), Seq(seq2))[0]
    dist = 0
    for i in range(len(alignment[0])):
        if alignment[0][i] != alignment[1][i]:
            dist += 1
    return dist

##########################################################################
def GetColorByNormalizedValue(cmap_name, norm_value):
    if norm_value < 0 or norm_value > 1:
        print("ERROR: value " + str(norm_value) + ' does not belong to [0, 1]')
    cmap = plt.cm.get_cmap(cmap_name)
    color = cmap(norm_value) 
    return mplt.colors.rgb2hex(color[:3])

##########################################################################
def GetBaseGeneName(gene_id):
    return gene_id.split('*')[0]

def PrepareOutputDir(output_dir):
    if os.path.exists(output_dir):
        shutil.rmtree(output_dir)
    os.mkdir(output_dir)

##########################################################################
def OutputPlotToPdf(output_fname):
    pp = PdfPages(output_fname)
    pp.savefig()
    pp.close()
    plt.clf()

##########################################################################
class StructuralRegion(Enum):
    FR1 = 0
    CDR1 = 1
    FR2 = 2
    CDR2 = 3
    FR3 = 4
    CDR3 = 5
    FR4 = 6
