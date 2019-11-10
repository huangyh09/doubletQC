## Command line for running scrublet
## https://pypi.org/project/scrublet

import os
import sys
import numpy as np
from scipy import io
import scrublet as scr
from optparse import OptionParser, OptionGroup


def load_10X(path, min_counts=None, min_cells=None, version3=False):
    """
    Load 10X data from cellranger output matrix, into 
    scipy csr matrix, arrays for genes and cell barcodes
    
    Filter cells by min_counts and filter genes by min_cells
    """
    ## load 10X matrix folder
    if version3:
        mat = io.mmread(path + "/matrix.mtx.gz").tocsr()
        genes = np.genfromtxt(path + "/features.tsv.gz", dtype="str", delimiter="\t")
        cells = np.genfromtxt(path + "/barcodes.tsv.gz", dtype="str", delimiter="\t")
    else:
        mat = io.mmread(path + "/matrix.mtx").tocsr()
        genes = np.genfromtxt(path + "/genes.tsv", dtype="str", delimiter="\t")
        cells = np.genfromtxt(path + "/barcodes.tsv", dtype="str", delimiter="\t")
    
    ## filter cells
    if min_counts is not None and min_counts > 0:
        n_counts = np.array(np.sum(mat, axis=0)).reshape(-1)
        idx = n_counts >= min_counts
        mat = mat[:, idx]
        cells = cells[idx]
       
    ## filter genes
    if min_cells is not None and min_cells > 0:
        n_cells = np.array(np.sum(mat, axis=1)).reshape(-1)
        idx = n_counts >= min_counts
        mat = mat[idx, :]
        genes = genes[idx, ]

    return mat, genes, cells


def main():
    # parse command line options
    parser = OptionParser()
    parser.add_option("--inputDir", "-i", dest="input_dir", default=None,
        help=("Directory of input matrix in 10x cellranger format"))
    parser.add_option("--outFile", "-o", dest="out_file", default=None,
        help=("Path for output file [default: $i/scrublet_table.tsv]"))
    parser.add_option("--cellranger2", "-2", dest="cellranger2", 
        action="store_true", default=False, 
        help="Use it for cellranger v2 instead of v3")
    parser.add_option("--expected_rate", "-r", dest="expected_rate", 
        default=None, help="Expected doublet rate: [default: n_cell/100K].")
    
    
    (options, args) = parser.parse_args()
    
    dat_path = os.path.abspath(options.input_dir)
    version3 = options.cellranger2 == False
    mat_dat, gene_ids, cell_ids = load_10X(dat_path, min_counts=None, 
                                           min_cells=None, version3=version3)
    
    n_cell = mat_dat.shape[1]
    if options.expected_rate is None:
        expected_rate = n_cell / 100000.0
    else:
        expected_rate = float(expected_rate)
    expected_rate = min(expected_rate, 0.5)
    
    print("Files loaded: %d cells." %(n_cell))
    print("Expected doublet rate: %.3f" %(expected_rate))
        
    scrub = scr.Scrublet(mat_dat.transpose(), 
                         expected_doublet_rate=expected_rate)
    scores, is_doublets = scrub.scrub_doublets()
    
    if options.out_file is None:
        out_file = dat_path + "/scrublet_table.tsv"
    else:
        out_file = options.out_file
    fid = open(out_file, "w")
    fid.writelines("cellID\tDBscore\tdoublet\n")
    for i in range(len(cell_ids)):
        out_list = [cell_ids[i], "%.3f" %scores[i], str(is_doublets[i])]
        fid.writelines("\t".join(out_list) + "\n")
    fid.close()
    
if __name__ == "__main__":
    main()
    