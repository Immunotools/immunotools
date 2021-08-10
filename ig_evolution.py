#!/usr/bin/env python2

import os
import shutil
import sys
import getopt

script_path = os.path.abspath(os.path.dirname(os.path.realpath(__file__)))
sys.path.append(os.path.join(script_path, 'py/ig_evolution'))
import utils
import dataset
import cdr3_clonal_lineage
import full_length_clonal_lineage
import clonal_tree_writer
import amino_acid_utils
import clonal_tree_utils
import clonal_tree_constructor
import clonal_graph_utils
import html_writer
import lineage_stats_writer

def FilterClonalLineage(cdr3_lineage, config):
    return cdr3_lineage.NumFullLengthSequences() < config.min_lineage_size or \
           cdr3_lineage.NumFullLengthSequences() > config.max_lineage_size

def PrepareOutputDirs(output_dir):
    utils.PrepareOutputDir(output_dir)
    output_dirs = dict()
    output_dirs['main_dir'] = output_dir
    # graph output
    output_dirs['clonal_graphs'] = os.path.join(output_dir, 'clonal_graphs')
    utils.PrepareOutputDir(output_dirs['clonal_graphs'])
    output_dirs['coloring_label'] = os.path.join(output_dirs['clonal_graphs'], 'labels')
    utils.PrepareOutputDir(output_dirs['coloring_label'])
    output_dirs['coloring_mult'] = os.path.join(output_dirs['clonal_graphs'], 'multiplicity')
    utils.PrepareOutputDir(output_dirs['coloring_mult'])    
    output_dirs['coloring_shm'] = os.path.join(output_dirs['clonal_graphs'], 'shm_depth')
    utils.PrepareOutputDir(output_dirs['coloring_shm'])
    output_dirs['compressed'] = os.path.join(output_dirs['clonal_graphs'], 'compressed')
    utils.PrepareOutputDir(output_dirs['compressed'])
    output_dirs['shm_matrix'] = os.path.join(output_dirs['clonal_graphs'], 'shm_matrix')
    utils.PrepareOutputDir(output_dirs['shm_matrix'])
    output_dirs['shm_plot'] = os.path.join(output_dirs['clonal_graphs'], 'shm_plot')
    utils.PrepareOutputDir(output_dirs['shm_plot'])
    # SHMs
    output_dirs['shms'] = os.path.join(output_dir, 'shms')
    utils.PrepareOutputDir(output_dirs['shms'])
    output_dirs['position_shms'] = os.path.join(output_dirs['shms'], 'shms_per_position')
    utils.PrepareOutputDir(output_dirs['position_shms'])
    # HTMLs
    output_dirs['htmls'] = os.path.join(output_dir, 'htmls')
    utils.PrepareOutputDir(output_dirs['htmls'])
    # aux dirs
    output_dirs['fl_lineages'] = os.path.join(output_dir, 'full_length_lineages')
    utils.PrepareOutputDir(output_dirs['fl_lineages'])
    output_dirs['cdr3_fasta_dir'] = os.path.join(output_dir, 'cdr3_fasta')
    utils.PrepareOutputDir(output_dirs['cdr3_fasta_dir'])
    output_dirs['cdr3_graph_dir'] = os.path.join(output_dir, 'cdr3_graphs')
    utils.PrepareOutputDir(output_dirs['cdr3_graph_dir'])
    output_dirs['cdr3_comp_dir'] = os.path.join(output_dir, 'cdr3_connected_components')
    utils.PrepareOutputDir(output_dirs['cdr3_comp_dir'])
    output_dirs['cdr3_lin_dir'] = os.path.join(output_dir, 'cdr3_lineages')
    utils.PrepareOutputDir(output_dirs['cdr3_lin_dir'])
    return output_dirs

def RemoveAuxDirs(output_dirs):
    shutil.rmtree(output_dirs['fl_lineages'])
    shutil.rmtree(output_dirs['cdr3_fasta_dir'])
    shutil.rmtree(output_dirs['cdr3_graph_dir'])
    shutil.rmtree(output_dirs['cdr3_comp_dir'])
    shutil.rmtree(output_dirs['cdr3_lin_dir'])

class AlgorithmConfig:
    def __init__(self, args):
        self.min_abs_abundance = 5
        self.max_rel_abundance = 0.01 #0.05
        self.min_graph_size = 10
        self.parse_headers = False
        self.min_lineage_size = 100
        self.max_lineage_size = sys.maxint
        self.hg_tau = 30
        self.min_component_fraction = 0.0
        self.perc_cdr3_identity = 90
        self.min_cdr3_len = 10
        self.divan_dir = ''
        self.output_dir = ''
        self.num_lineages = sys.maxint
        self.remove_aux_files = True
        self._ParseArgs(args)
        self._PrintArgs()

    def _PrintArgs(self):
        print "==== IgEvolution parameters:"
        print "Min relative (absolute) vertex abundance: " + str(self.max_rel_abundance) + ' (' + str(self.min_abs_abundance) + ')'
        print "Min lineage (graph) size: " + str(self.min_lineage_size) + ' (' + str(self.min_graph_size) + ')'
        print "Hamming graph TAU: " + str(self.hg_tau) 

    def _ParseArgs(self, args):
        try:
            options, remainder = getopt.getopt(args[1:], 'i:o:n:', ["parse-mults", 'min-lineage=', 'max-lineage=', 'min-abs=',
                                                                  'min-rel=', 'hg-tau=', 'min-graph=', 'help', 'num-lineages=',
                                                                  'skip-err-corr', 'keep-aux-files', 'all', 'cdr3-pi=', 'min-comp-fr=', 'min-cdr3='])
        except getopt.GetoptError as err:
            print str(err)  # will print something like "option -a not recognized"
            sys.exit(2)
        for opt, arg in options:
            if opt == "-i":
                self.divan_dir = arg
            elif opt == '-o':
                self.output_dir = arg
            elif opt == '--parse-mults':
                self.parse_headers = True
            elif opt == '--min-lineage':
                self.min_lineage_size = int(arg)
            elif opt == '--max-lineage':
                self.max_lineage_size = int(arg)
            elif opt == '--min-abs':
                self.min_abs_abundance = int(arg)
            elif opt == '--min-rel':
                self.max_rel_abundance = float(arg)
            elif opt == '--min-graph':
                self.min_graph_size = int(arg)
            elif opt == '--hg-tau':
                self.hg_tau = int(arg)
            elif opt == '--skip-err-corr':
                self.min_abs_abundance = 0
                self.max_rel_abundance = 0
            elif opt == '--keep-aux-files':
                self.remove_aux_files = False
            elif opt == '--cdr3-pi':
                self.perc_cdr3_identity = float(arg)
            elif opt == '--min-comp-fr':
                self.min_component_fraction = float(arg)
            elif opt == '--min-cdr3':
                self.min_cdr3_len = int(arg)
            elif opt == '--all':
                self.min_graph_size = 1
                self.min_lineage_size = 1
            elif opt == '-n' or opt == '--num-lineages':
                self.num_lineages = int(arg)
            elif opt == '--help':
                self._PrintHelp()
                sys.exit(0)
            else:
                assert False, "unhandled option"
        if self.output_dir == '':
            print "ERROR: output directory (-o) was not specified"
            sys.exit(1)
        if self.divan_dir == '':
            print "ERROR: input directory (-i) was not specified"
            sys.exit(1)

    def _PrintHelp(self):
        print "python ig_evolution_launch.py -i divan_dir -o output_dir [--min-lineage MIN_LIN_SIZE " \
              "--max-lineage MAX_LIN_SIZE --min-graph MIN_GRAPH_SIZE --hg-tau HG_TAU --min-abs MIN_ABS_ABUN " \
              "--min-rel MIN_REL_ABUN --parse-mults --skip-err-corr]"

def main(argv):
    config = AlgorithmConfig(argv)
    divan_test_dir = config.divan_dir
    output_dir = config.output_dir

    output_dirs = PrepareOutputDirs(output_dir)

    dataset_info = dataset.DatasetInfo('test1', 'time_point1', 1, 'cell_type1')
    dataset_obj = dataset.Dataset(dataset_info, divan_test_dir, config.parse_headers)

    cdr3_lineage_constructor = cdr3_clonal_lineage.CDR3LineageConstructor(dataset_obj, output_dirs, config.perc_cdr3_identity)
    cdr3_lineages = cdr3_lineage_constructor.Construct()

    full_length_lineages = []
    for l in cdr3_lineages:
        full_length_lineages.append(full_length_clonal_lineage.FullLengthClonalLineage(l))
    print str(len(full_length_lineages)) + " full-length lineages were constructed"

    print "Writing clonal lineage statistics..."
    lineage_stats = lineage_stats_writer.ClonalLineageStatWriter(full_length_lineages)
    lineage_stats.OutputStats(os.path.join(output_dirs['main_dir'], 'raw_lineage_stats.txt'))

#    clonal_graph_utils.OutputAbundantAAGraphs(full_length_lineages, output_dirs, config)
#
#    print "Compiling HTML report..."
#    dir_dict = {'labels' : output_dirs['coloring_label'], 'multiplicity' : output_dirs['coloring_mult'], 'compressed' : output_dirs['compressed'], 'shm_matrix' : output_dirs['shm_matrix'], 'shm_plot' : output_dirs['shm_plot'], 'shm_depth' : output_dirs['coloring_shm']}
#    html = html_writer.HTMLWriter(os.path.basename(output_dirs['clonal_graphs']), dir_dict, '.svg', ['shm_matrix', 'shm_plot'])
#    html.CreateHTMLReports(output_dirs['htmls'])

    if config.remove_aux_files:
        print "Removing auxiliary direcitories..."
        RemoveAuxDirs(output_dirs)

    print "Thank you for using IgEvolution!"

if __name__ == '__main__':
    main(sys.argv)
