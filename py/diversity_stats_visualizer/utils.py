import logging
import sys
import imp

import matplotlib as mplt
mplt.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

def get_logger_by_arg(logger_arg, logger_name = ""):
    if logger_arg != "-":
        return logger_arg
    log = logging.getLogger(logger_name)
    log.setLevel(logging.DEBUG)
    console = logging.StreamHandler(sys.stdout)
    console.setFormatter(logging.Formatter('%(message)s'))
    console.setLevel(logging.DEBUG)
    log.addHandler(console)
    return log

def output_figure(output_fname, log_info, log):
    pp = PdfPages(output_fname + ".pdf")
    pp.savefig()
    pp.close()
    plt.savefig(output_fname + ".png")
    plt.savefig(output_fname + ".svg")
    plt.clf()
    log.info(log_info + " was written to " + output_fname + ".pdf, .png and .svg")

hydrophobic = ['L', 'I', 'F', 'W', 'V', 'M', 'C', 'Y', 'A']
neutral = ['T', 'E', 'G', 'S', 'Q', 'D', '*']
hydrophilic = ['R', 'K', 'N', 'H', 'P']

hydrophoby_dict = {'L' : 100, 'I' : 100, 'F' : 92, 'W' : 84, 'V' : 79, 'M' : 74, 'C' : 52, 'Y' : 49, 'A' : 47,
                   'T' : 13, 'E' : 8, 'G' : 0, 'S' : -7, 'Q' : -18, 'D' : -18, '*' : 0, 'X' : 0,
                   'R' : -26, 'K' : -37, 'N' : -41, 'H' : -42, 'P' : -46}

hydro_min = -100
hydro_max = 100

###################################################
def CheckPackageFatal(package_name, log):
    try:
        imp.find_module(package_name)
        found = True
    except ImportError:
        found = False
    if not found:
        log.info("ERROR: package " + package_name + ' is not installed. Please look manual for installation details')
        sys.exit(1)

def GetBaseName(gene_name):
    return gene_name.split('*')[0] 
