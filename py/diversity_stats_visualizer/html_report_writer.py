import os
import sys
import warnings

import pandas as pd

import visualize_shm_stats
import visualize_cdr_stats
import utils

class HTMLReportWriter:
    def _WriteHeader(self, level, text, align):
        self.fhandler.write("<h" + str(level) + " align = " + align + ">" + text + "</h" + str(level) + ">\n")

    def _GetAlignment(self, align_to_center = True):
        if not align_to_center:
            return "left"
        return "center"

    def _WriteTableCell(self, elem):
        self.fhandler.write("<td>" + str(elem) + "</td>")

    def _WriteTableRow(self, row):
        self.fhandler.write("<tr>\n")
        for r in row:
            self._WriteTableCell(r)
        self.fhandler.write("</tr>\n")

    def __init__(self, output_fname):
        self.fhandler = open(output_fname, "w")

    def WriteH1(self, text, align_to_center = True):
        self._WriteHeader(1, text, self._GetAlignment(align_to_center))

    def WriteH2(self, text, align_to_center = True):
        self._WriteHeader(2, text, self._GetAlignment(align_to_center))

    def WriteH3(self, text, align_to_center = True):
        self._WriteHeader(3, text, self._GetAlignment(align_to_center))

    # width in percent
    def WriteImage(self, path_to_image, width = 60):
        self.fhandler.write("<p align = center>\n")
        self.fhandler.write("<image src = " + path_to_image + " width = " + str(width) + "%></image>\n")
        self.fhandler.write("</p>\n")

    def WriteTable(self, col_names, row_names, values, width = 60):
        if len(values) == 0:
            return
        if len(row_names) != len(values):
            print "# rows in table and # row names are not consistent"
            sys.exit(1)
        if len(col_names) != len(values[0]):
            print "# columns in table and # column names are not consistent"
            sys.exit(1)
        self.fhandler.write("<table width = " + str(width) + "% align = center>\n")
        header_row = [""]
        header_row.extend(col_names)
        self._WriteTableRow(header_row)
        for i in range(0, len(row_names)):
            cur_row = [row_names[i]]
            cur_row.extend(values[i])
            self._WriteTableRow(cur_row)
        self.fhandler.write("</table>\n")

    def WriteEmptyLine(self):
        self.fhandler.write("<br>\n")

    def WriteHorizontalLine(self, width = 100):
        self.fhandler.write("<hr width=" + str(width) + "%>\n")

    def CloseFile(self):
        self.fhandler.close()

    def WriteImageWithTitle(self, fname, title, output_config):
        if os.path.exists(fname):
            self.WriteH2(title)
            self.WriteImage(output_config.GetFnameForHTML(fname), width=60)    

#######################################################################################################################
def ComputeLocusCharacteristics(vj_df, locus):
    stats = []
    locus_df = vj_df.loc[vj_df['Chain_type'] == locus]
    stats.append(len(locus_df["Read_name"]))
    stats.append(len([s for s in locus_df['Productive'] if s == 1]))
    stats.append(len([s for s in locus_df['Has_stop_codon'] if s == 1]))
    stats.append(len([s for s in locus_df['In-frame'] if s == 0]))
    return stats

def ComputeGeneralCharacteristicsTable(vj_df):
    col_names = ["IGH", "IGK", "IGL"]
    row_names = ["# sequences", "# productive sequences", "# sequences with stop codon", "# out-of-frame sequences"]
    table = []
    for r in row_names:
        table.append([0] * len(col_names))
    for i in range(0, len(col_names)):
        locus_stats = ComputeLocusCharacteristics(vj_df, col_names[i])
        for j in range(0, len(locus_stats)):
            table[j][i] = locus_stats[j]
    return table, row_names, col_names

def WriteGeneralCharacteristics(html_writer, vj_df, output_config):
    html_writer.WriteH1("General characteristics")
    table, row_names, col_names = ComputeGeneralCharacteristicsTable(vj_df)
    html_writer.WriteTable(col_names, row_names, table)
    html_writer.WriteImageWithTitle(output_config.vj_usage + '.svg', "Heatmap of VJ hit concentrations", output_config)
    html_writer.WriteImageWithTitle(output_config.v_usage + '.svg', "V usages", output_config)
    html_writer.WriteImageWithTitle(output_config.j_usage + '.svg', "J usages", output_config)

#######################################################################################################################
def ComputeGeneralSHMCharacteristics(shm_df):
    col_names = ['V segment', 'J segment']
    row_names = ['# SHMs', '# substitutions', '# deletions', '# insertions', '# synonymous SHMs', '# stop codon SHMs']
    row_index_dict = dict()
    for i in range(0, len(row_names)):
        row_index_dict[row_names[i]] = i
    table = []
    for r in row_names:
        table.append([0] * len(col_names))
    for read in shm_df:
        shms = shm_df[read]
        index_j = 0
        if not read.is_variable():
            index_j = 1
        for shm in shms:
            table[row_index_dict["# SHMs"]][index_j] += 1
            if shm.is_substitution():
                table[row_index_dict["# substitutions"]][index_j] += 1
            if shm.is_deletion():
                table[row_index_dict["# deletions"]][index_j] += 1
            if shm.is_insertion():
                table[row_index_dict["# insertions"]][index_j] += 1
            if shm.synonymous:
                table[row_index_dict["# synonymous SHMs"]][index_j] += 1
            if shm.to_stop_codon:
                table[row_index_dict["# stop codon SHMs"]][index_j] += 1
    return row_names, col_names, table

def WriteSHMCharacteristics(html_writer, shm_df, output_config):
    html_writer.WriteH1("SHM characteristics")
    row_names, col_names, table = ComputeGeneralSHMCharacteristics(shm_df)
    html_writer.WriteTable(col_names, row_names, table)
    html_writer.WriteImageWithTitle(output_config.aa_matrix + '.svg', "Heatmap of amino acid substitutions", output_config)
    html_writer.WriteImageWithTitle(output_config.nucl_matrix + '.svg', "Heatmap of nucleotide substitutions", output_config)
    for l, fname in output_config.NumSHMIter():
        html_writer.WriteImageWithTitle(fname + '.svg', 'Distribution of SHM in ' + l + 'V', output_config)
    html_writer.WriteImageWithTitle(output_config.indel_length + '.svg', 'Distribution of lengths of insertion/deletion V SHMs', output_config)

#######################################################################################################################
def ComputeLocusCDRCharacteristics(vj_df, locus):
    locus_df = vj_df.loc[vj_df['Chain_type'] == locus]
    nums = list()
    cdrs = ['CDR1', 'CDR2', 'CDR3']
    for cdr in cdrs:
        region_seq = list(locus_df[cdr + "_nucls"])
        seq_dict = dict()
        for seq in region_seq:
            if seq not in seq_dict:
                seq_dict[seq] = 0
            seq_dict[seq] += 1
        max_abun = 0
        for seq in seq_dict:
            max_abun = max(max_abun, seq_dict[seq])
        nums.extend([len(seq_dict), max_abun])
    return nums

def ComputeGeneralCDRCharacteristicsTable(vj_df):
    col_names = ["IGH", "IGK", "IGL"]
    row_names = ["# CDR1s", "max CDR1 abundance",
                 "# CDR2s", "max CDR2 abundance",
                 "# CDR3s", "max CDR3 abundance"]
    table = []
    for r in row_names:
        table.append([0] * len(col_names))
    for i in range(0, len(col_names)):
        nums = ComputeLocusCDRCharacteristics(vj_df, col_names[i])
        for j in range(0, len(nums)):
            table[j][i] = nums[j]
    return table, row_names, col_names

def WriteGeneralCDRCharacteristics(html_writer, vj_df):
    table, row_names, col_names = ComputeGeneralCDRCharacteristicsTable(vj_df)
    html_writer.WriteTable(col_names, row_names, table)

def GetLengthAbundanceLargestGroup(locus_df, locus, cdr):
    region_seq = list(locus_df[cdr + "_nucls"])
    max_group = visualize_cdr_stats.get_region_largest_group(region_seq)
    if len(max_group) == 0:
        return 0, 0
    return len(max_group[0]), round(float(len(max_group)) / float(len(region_seq)) * 100, 2)

def ImageDictContainsLocus(image_dict, locus):
    for im in image_dict:
        if im[:len(locus)] == locus:
            return True
    return False

def WriteCDRPlots(html_writer, vj_df, images_dict):
    loci = ['IGH', 'IGK', 'IGL']
    cdrs = ['CDR1', 'CDR2', 'CDR3']
    for l in loci:
        if not ImageDictContainsLocus(images_dict, l):
            continue
        locus_df = vj_df.loc[vj_df['Chain_type'] == l]
        for cdr in cdrs:
            if l + "_" + cdr + "_length" in images_dict:
                html_writer.WriteH2(l + " " + cdr + " length distribution:")
                html_writer.WriteImage(images_dict[l + "_" + cdr + "_length"], 65)
            else:
                continue
            max_len, max_abun = GetLengthAbundanceLargestGroup(locus_df, l, cdr)
            if l + "_" + cdr + "_nucls" in images_dict:
                html_writer.WriteH2("Most " + l + " " + cdr + "s (" + str(max_abun) + "%) have length " + str(max_len) + " nt")
                html_writer.WriteH3("Distribution of nucleotide abundance per position for " + l + " " + cdr +
                                    "s of length " + str(max_len) + " nt:")
                html_writer.WriteImage(images_dict[l + "_" + cdr + "_nucls"], 70)
                if l + "_" + cdr + "_aa" in images_dict:
                    html_writer.WriteH3(l + " " + cdr + "s of length " + str(max_len) + " nt are in-frame. "
                                                                                   "Plot of the most abundant amino "
                                                                                   "acids per position:")
                    html_writer.WriteImage(images_dict[l + "_" + cdr + "_aa"], 70)
        html_writer.WriteHorizontalLine(70)

def WriteCDRCharacteristics(html_writer, vj_df, images_dict):
    html_writer.WriteH1("CDR characteristics")
    WriteGeneralCDRCharacteristics(html_writer, vj_df)
    WriteCDRPlots(html_writer, vj_df, images_dict)

#######################################################################################################################
def create_html(vj_df, shm_df, output_config):
    output_config.Log().info("Annotation report will be written to " + output_config.html_report)
    output_config.Log().info("Printing general characteristics of the repertoire")
    html_writer = HTMLReportWriter(output_config.html_report)
    WriteGeneralCharacteristics(html_writer, vj_df, output_config)
    html_writer.WriteHorizontalLine()
    output_config.Log().info("Printing SHM characteristics")
    WriteSHMCharacteristics(html_writer, shm_df, output_config)
    html_writer.WriteHorizontalLine()
    output_config.Log().info("Printing CDR characteristics")
#    WriteCDRCharacteristics(html_writer, vj_df, images_dict)
    html_writer.CloseFile()
    output_config.Log().info("Annotation report was written to " + output_config.html_report)

#######################################################################################################################
def main(vj_df_fname, shm_df_fname, output_config):
    if not os.path.exists(vj_df_fname):
        output_config.Log().info("CDR details file " + vj_df_fname + " was not found")
        sys.exit(1)
    vj_df = pd.read_table(vj_df_fname, delim_whitespace = True)
    if not os.path.exists(shm_df_fname):
        output_config().Log().info("SHM details file " + shm_df_fname + " was not found")
        sys.exit(1)
    shm_df = visualize_shm_stats.SHMs(shm_df_fname)
    create_html(vj_df, shm_df, output_config)

if __name__ == "__main__":
    warnings.filterwarnings('ignore')
    if len(sys.argv) != 4:
        print "Invalid input"
        print "python html_report_writer.py cdr_details.txt shm_details.txt output_config"
        sys.exit(1)
    main(sys.argv[1], sys.argv[2], sys.argv[3])
