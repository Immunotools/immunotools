import os
import sys

import utils

import pylab
import matplotlib as mplt
mplt.use('Agg')

import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.patches as patches

class AlignedRead:
    def _parse_line(self, line):
        splits = line.split()
        if len(splits) != 6:
            print "Header line " + line + " is not valid"
            sys.exit(1)
        self.read_name = splits[0][len("Read_name:"):]
        self.read_len = int(splits[1][len("Read_length:"):])
        self.gene_name = splits[2][len("Gene_name:"):]
        self.gene_len = int(splits[3][len("Gene_length:"):])
        self.segment = splits[4][len("Segment:"):]
        self.chain_type = splits[5][len("Chain_type:"):]

    def __init__(self, line):
        self._parse_line(line)

    def __hash__(self):
        return hash(self.read_name) * hash(self.gene_name)

    def __eq__(self, other):
        return self.read_name == other.read_name and self.gene_name == other.gene_name

    def is_variable(self):
        return self.segment == 'V'

class SHM:
    def __init__(self, line):
        splits = line.split()
        if len(splits) != 9:
            print "Invalid SHM line " + line
        self.type = splits[0]
        self.read_pos = int(splits[1])
        self.gene_pos = int(splits[2])
        self.read_nucl = splits[3]
        self.gene_nucl = splits[4]
        self.read_aa = splits[5]
        self.gene_aa = splits[6]
        self.synonymous = False
        if int(splits[7]) == 1:
            self.synonymous = True
        self.to_stop_codon = False
        if int(splits[8]) == 1:
            self.to_stop_codon = True

    def is_deletion(self):
        return self.type == 'D'

    def is_insertion(self):
        return self.type == 'I'

    def is_substitution(self):
        return self.type == 'S'

class SHMs:
    def __init__(self, df_fname):
        self.shm_dict = dict()
        fhandler = open(df_fname, "r")
        lines = fhandler.readlines()
        current_read = ""
        for i in range(1, len(lines)):
            l = lines[i].strip()
            if l[:len("Read_name:")] == "Read_name:":
                current_read = AlignedRead(l)
                self.shm_dict[current_read] = list()
            else:
                self.shm_dict[current_read].append(SHM(l))

    def __len__(self):
        return len(self.shm_dict)

    def __getitem__(self, item):
        return self.shm_dict[item]

    def __iter__(self):
        for it in self.shm_dict:
            yield it

isotype_colors = {'IGH': 'b', 'IGK': 'g', 'IGL': 'r'}

############################### NUMBER of SHMs per GENE #################################################
def OutputGeneSHMPlot(gene_shms, gene_name, gene_length, num_aligned_seqs, output_fname, log):
    nucl_dict = {'A' : [0] * gene_length, 'C' : [0] * gene_length, 'G' : [0] * gene_length, 'T' : [0] * gene_length}
    num_shms = 0
    for shm in gene_shms:
        if not shm.is_substitution() or not nucl_is_valid(shm.read_nucl):
            continue
        nucl_dict[shm.read_nucl][shm.gene_pos] += 1
        num_shms += 1
    for nucl in nucl_dict:
        for i in range(len(nucl_dict[nucl])):
            nucl_dict[nucl][i] = nucl_dict[nucl][i]
    x = range(gene_length)
    plt.figure()
    plt.bar(x, [float(sum(y)) / num_aligned_seqs for y in zip(nucl_dict['A'], nucl_dict['C'], nucl_dict['G'], nucl_dict['T'])], color = 'blue', label = 'A')
    plt.bar(x, [float(sum(y)) / num_aligned_seqs for y in zip(nucl_dict['A'], nucl_dict['C'], nucl_dict['G'])], color = 'green', label = 'C')
    plt.bar(x, [float(sum(y)) / num_aligned_seqs for y in zip(nucl_dict['A'], nucl_dict['C'])], color = 'red', label = 'G')
    plt.bar(x, [float(m) / num_aligned_seqs for m in nucl_dict['T']], color = 'orange', label = 'T')
    plt.legend(loc = 'upper center', ncol = 4)
    plt.ylim(0, 1.1)
    plt.xlabel('Position in V gene', fontsize = 14)
    plt.ylabel('Fraction of sequences', fontsize = 14)
    plt.title(str(num_aligned_seqs) + ' sequences were aligned to ' + gene_name)
    utils.output_figure(output_fname, "SHM position in " + gene_name, log)
    return nucl_dict

def WriteNucleotide(fh, pos, nucl, nucl_mult, num_aligned_seq):
    if nucl_mult == 0:
        return
    fh.write(str(pos) + '\t' + nucl + '\t' + str(nucl_mult) + '\t' + str(float(nucl_mult) / num_aligned_seq) + '\n')

def OutputGeneSHMsToTxt(nucl_dict, num_aligned_seqs, output_fname):
    fh = open(output_fname, 'w')
    fh.write('Position\tNucleotide\tMultiplicity\tFrequency\n')
    for i in range(len(nucl_dict['A'])):
        for nucl in ['A', 'C', 'G', 'T']:
            WriteNucleotide(fh, i, nucl, nucl_dict[nucl][i], num_aligned_seqs)
    fh.close()

def OutputSHMsForVGenes(shm_df, output_config):
    gene_type_dict = dict()
    gene_len = dict()
    num_aligned = dict()
    for it in shm_df:
        if it.segment not in gene_type_dict:
            gene_type_dict[it.segment] = dict()
        segment_dict = gene_type_dict[it.segment]
        gene_name = utils.GetBaseName(it.gene_name)
        if gene_name not in gene_type_dict[it.segment]:
            gene_type_dict[it.segment][gene_name] = []
            gene_len[gene_name] = 0
            num_aligned[gene_name] = 0
        gene_type_dict[it.segment][gene_name].extend(shm_df[it])
        gene_len[gene_name] = max(gene_len[gene_name], it.gene_len)
        num_aligned[gene_name] += 1
    for segment in gene_type_dict:
        segment_dict = gene_type_dict[segment]
        for gene_name in segment_dict:
            num_aligned_seq = num_aligned[gene_name]
            if num_aligned_seq < 10:
                continue
            output_fname = os.path.join(output_config.GetSHMDirBySegment(segment), gene_name)
            nucl_pos_dict = OutputGeneSHMPlot(segment_dict[gene_name], gene_name, gene_len[gene_name], num_aligned[gene_name], output_fname, output_config.Log())
            output_config.AddSHMFileForSegment(segment, output_fname)
            OutputGeneSHMsToTxt(nucl_pos_dict, num_aligned[gene_name], os.path.join(output_config.GetSHMDirBySegment(segment), gene_name) + '.txt')

def OutputGeneMutability(gene_mutability_dict, output_fname, gene_type, log):
    df_dict = {'Gene' : [], 'Mutability' : []}
    for gene in gene_mutability_dict:
        for m in gene_mutability_dict[gene]:
            df_dict['Gene'].append(gene)
            df_dict['Mutability'].append(m)
    plt.figure(figsize = (10, 8))
    sns.boxplot(x = 'Gene', y = 'Mutability', data = df_dict)
    max_mutability = max(0.55, max(df_dict['Mutability']))
    plt.ylim(-0.05, max_mutability)
    plt.xticks(rotation = 90)
    plt.ylabel('Mutability')
    utils.output_figure(output_fname, "Mutability of " + gene_type + ' genes', log)

def OutputVJGenesMutability(shm_df, output_config):
    v_gene_mutability = dict()
    j_gene_mutability = dict()
    for it in shm_df:
        cur_dict = v_gene_mutability
        if not it.is_variable():
            cur_dict = j_gene_mutability
        gene_name = utils.GetBaseName(it.gene_name)
        if gene_name not in cur_dict:
            cur_dict[gene_name] = []
        mutability = float(len(shm_df[it])) / it.gene_len
        cur_dict[gene_name].append(mutability)
    OutputGeneMutability(v_gene_mutability, output_config.v_mutability, 'V', output_config.Log())
    OutputGeneMutability(j_gene_mutability, output_config.j_mutability, 'J', output_config.Log())

############################### NUMBER of SHMs per ISOTYPE #################################################
def output_shm_stats_for_isotype(num_shms, locus, output_fname, log):
    plt.figure()
    plt.hist(num_shms, color = isotype_colors[locus], bins = 50, alpha = .75)
    plt.xlabel("# SHMs in " + locus + "V", fontsize = 16)
    plt.ylabel("# sequences", fontsize = 16)
    plt.xticks(fontsize = 14)
    plt.yticks(fontsize = 14)
    plt.title('# SHMs in ' + locus + 'V' + ' sequences', fontsize = 14)
    utils.output_figure(output_fname, "Distribution of # SHMs in " + locus + "V segments", log)

def ComputeNumSHMsInLoci(shm_df):
    locus_dict = dict()
    for it in shm_df:
        if not it.is_variable():
            continue
        if it.chain_type not in locus_dict:
            locus_dict[it.chain_type] = []
        locus_dict[it.chain_type].append(len(shm_df[it]))
    return locus_dict

def visualize_v_mutations_stats(shms_df, output_config):
    num_shms_in_loci = ComputeNumSHMsInLoci(shms_df)
    for l, fname in output_config.NumSHMIter():
        num_shms = 0
        if l in num_shms_in_loci:
            num_shms = len(num_shms_in_loci[l])
        if num_shms < 10:
            output_config.Log().info("# sequences for " + l + " is too small (" + str(num_shms) + "). Plot drawing was skipped")
            continue
        output_shm_stats_for_isotype(num_shms_in_loci[l], l, fname, output_config.Log())

#################################### SUBSTITUTION MATRICES ###################################################
def get_aa_list():
    return ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']

def aa_is_valid(aa):
    return aa != '*' and aa != '-' and aa != 'X'

def get_aa_ticks_colors(aa_list):
    colors = []
    for aa in aa_list:
        if aa in utils.hydrophobic:
            colors.append('red')
        elif aa in utils.neutral:
            colors.append('green')
        elif aa in utils.hydrophilic:
            colors.append('blue')
    return colors

def visualize_aa_substitution_matrix(shms_df, output_fname, log):
    dict_aa = dict()
    num_shms = 0
    for it in shms_df:
        read_shms = shms_df[it]
        prev_pos = -1
        for shm in read_shms:
            if prev_pos / 3 != shm.read_pos / 3:
                if aa_is_valid(shm.gene_aa) and aa_is_valid(shm.read_aa):
                    aa_pair = shm.gene_aa + shm.read_aa
                    if not aa_pair in dict_aa:
                        dict_aa[aa_pair] = 0
                    dict_aa[aa_pair] += 1
                    num_shms += 1
            prev_pos = shm.read_pos
    aa_list = get_aa_list()
    aa_freq = []
    for i in range(0, len(aa_list)):
        aa_freq.append([0] * len(aa_list))
    for aa_pair in dict_aa:
        aa_freq[aa_list.index(aa_pair[1])][aa_list.index(aa_pair[0])] = float(dict_aa[aa_pair]) / float(num_shms)
    fig, ax = plt.subplots()
    sns.heatmap(aa_freq, cmap = plt.cm.jet, xticklabels = aa_list, yticklabels = aa_list, square = True, ax = ax)
    ax.tick_params(labelsize = 14)
    plt.xticks(fontsize = 12)
    plt.yticks(fontsize = 12, rotation='horizontal')
    #for tick, color in zip(ax.get_xticklabels(), get_aa_ticks_colors(aa_list)):
    #    tick.set_color(color)
    #for tick, color in zip(ax.get_yticklabels(), get_aa_ticks_colors(aa_list)):
    #    tick.set_color(color)
    plt.xlabel("To", fontsize = 14)
    plt.ylabel("From", fontsize = 14, rotation='horizontal')
    utils.output_figure(output_fname, "Amino acid substitution heatmap", log)
    return aa_freq

def nucl_is_valid(nucl):
    return nucl != 'N'

def visualize_nucl_substitution_matrix(shms_df, output_fname, log):
    nucl_list = ['A', 'C', 'G', 'T']
    nucl_matrix = []
    for n in nucl_list:
        nucl_matrix.append([0] * len(nucl_list))
    num_shms = 0
    for it in shms_df:
        read_shms = shms_df[it]
        for shm in read_shms:
            if not shm.is_substitution():
                continue
            if nucl_is_valid(shm.read_nucl) and nucl_is_valid(shm.gene_nucl):
                nucl_matrix[nucl_list.index(shm.read_nucl)][nucl_list.index(shm.gene_nucl)] += 1
                num_shms += 1
    for i in range(0, len(nucl_matrix)):
        for j in range(0, len(nucl_matrix[i])):
            nucl_matrix[i][j] = float(nucl_matrix[i][j]) / float(num_shms)
    fig, ax = plt.subplots()
    sns.heatmap(nucl_matrix, cmap = plt.cm.Blues, xticklabels = nucl_list, yticklabels = nucl_list, square = True, ax = ax)
    ax.tick_params(labelsize = 14)
    plt.xticks(fontsize = 12)
    plt.yticks(fontsize = 12, rotation='horizontal')
    plt.xlabel("To", fontsize = 14)
    plt.ylabel("From", fontsize = 14, rotation='horizontal')
    utils.output_figure(output_fname, "Nucleotide substitution heatmap", log)

def output_aa_freq(aa_freq, output_fname, log):
    aa_list = get_aa_list()
    fhandler = open(output_fname, "w")
    fhandler.write("from/to\t" + "\t".join(aa_list) + "\n")
    for i in range(0, len(aa_list)):
        fhandler.write(aa_list[i] + "\t" + "\t".join([str(ff) for ff in aa_freq[i]]) + "\n")
    log.info("Amino acid substitution matrix was written to " + output_fname)

#################################### SPECIAL SHMs ###################################################
def output_synonymous_shms(synonymous_pos, output_fname, log):
    if len(synonymous_pos) < 100:
        return
    plt.hist(synonymous_pos, color = 'r', bins = 100)
    plt.xlabel("Relative position of V SHM in read", fontsize = 14)
    plt.ylabel("#SHMs", fontsize = 14)
    plt.xlim(0, .75)
    plt.xticks(fontsize = 12)
    plt.yticks(fontsize = 12)
    utils.output_figure(output_fname, "Distribution of synonymous SHM positions in V segment", log)

def visualize_indel_shm_lengths(shm_df, output_fname, log):
    prev_read_pos = -1
    prev_gene_pos = -1
    insertion_length = []
    deletions_lengths = []
    in_len = 0
    del_len = 0
    for it in shm_df:
        read_shms = shm_df[it]
        for shm in read_shms:
            if shm.is_deletion():
                if shm.gene_pos - prev_gene_pos == 1:
                    del_len += 1
                else:
                    if del_len > 0:
                        deletions_lengths.append(del_len)
                    del_len = 1
                prev_gene_pos = shm.gene_pos
            if shm.is_insertion():
                if shm.read_pos - prev_read_pos == 1:
                    in_len += 1
                else:
                    if in_len > 0:
                        insertion_length.append(in_len)
                    in_len = 1
                prev_read_pos = shm.read_pos
    if in_len != 0:
        insertion_length.append(in_len)
    if del_len != 0:
        deletions_lengths.append(del_len)
    dt = []
    labels = []
    max_x_value = 0
    if len(deletions_lengths) > 10:
        dt.append(deletions_lengths)
        labels.append("Deletions")
    if len(insertion_length) > 10:
        dt.append(insertion_length)
        labels.append("Insertions")
    if len(dt) == 0:
        log.info("Output contains very low number of indel SHMs. Plot drawing was skipped")
        return
    plt.hist(dt, label = labels, bins = 50)
    plt.legend(loc = 'upper center', ncol = len(dt), fontsize = 14)
    plt.xlabel("Insertion / deletion SHM length", fontsize = 16)
    plt.ylabel("# insertion / deletion SHMs", fontsize = 16)
    xlim_right = 0
    if len(deletions_lengths) != 0:
        xlim_right = max(deletions_lengths)
    if len(insertion_length) != 0:
        xlim_right = max(xlim_right, max(insertion_length)) 
    plt.xlim(.5, xlim_right + .5)
    plt.xticks(range(0, xlim_right + 1), fontsize = 14)
    plt.yticks(fontsize = 14)
    utils.output_figure(output_fname, "Distribution of insertion/deletion SHM lengths", log)

def OutputFractionOfSynonymousSHMs(shm_df, output_fname, log):
    v_shm_fractions = []
    j_shm_fractions = []
    for it in shm_df:
        read_shms = shm_df[it]
        num_synonymous = 0
        for shm in read_shms:
            if shm.synonymous:
                num_synonymous += 1
        fraction = 0
        if len(read_shms) != 0:
            fraction = float(num_synonymous) / len(read_shms)
        if it.is_variable():
            v_shm_fractions.append(fraction)
        else:
            j_shm_fractions.append(fraction)
    plt.hist([v_shm_fractions, j_shm_fractions], label = ['V gene', 'J gene'])
    plt.xlabel('Fraction of synonymous SHMs', fontsize = 14)
    plt.ylabel('# sequences', fontsize = 14)
    plt.legend(loc = 'upper right', fontsize = 14)
    utils.output_figure(output_fname, "Fractions of synonymous SHMs in V and J genes", log)

#################################### MAIN ###################################################
def main(shm_df_fname, output_config):
    shm_df = SHMs(shm_df_fname)
    output_config.Log().info(str(len(shm_df)) + " records were extracted from " + shm_df_fname)
    if len(shm_df) == 0:
        output_config.Log().info("SHM data-frame contains 0 records. SHM visualization will be skipped")
        return
    visualize_v_mutations_stats(shm_df, output_config)
    OutputSHMsForVGenes(shm_df, output_config)
    OutputVJGenesMutability(shm_df, output_config)
    aa_freq = visualize_aa_substitution_matrix(shm_df, output_config.aa_matrix, output_config.Log())
#    output_aa_freq(aa_freq, os.path.join(output_dir, "aa_substitution_matrix.txt"), log)
    visualize_nucl_substitution_matrix(shm_df, output_config.nucl_matrix, output_config.Log())
    visualize_indel_shm_lengths(shm_df, output_config.indel_length, output_config.Log())
    OutputFractionOfSynonymousSHMs(shm_df, output_config.synonymous_shms, output_config.Log())

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print "Invalid input parameters"
        print "python visualize_shm_stats.py shm_df.txt output_config"
        sys.exit(1)
    main(sys.argv[1], sys.argv[2])
