## Command line for combining results of scrublet and doubletFinder
## https://github.com/huangyh09/doubletQC

## Author: Yuanhua Huang
## Date: November, 2019

import os
import sys
import numpy as np
import pandas as pd
from optparse import OptionParser, OptionGroup


def average_score(x1, x2, x3, W=np.array([0.25, 0.25, 0.5])):
    """Average scores with weights"""
    X = np.zeros((len(x1), 3))
    X[:, 0] = x1
    X[:, 1] = x2
    X[:, 2] = x3
    
    if W is None:
        W = np.ones((X.shape[1], 1)) / X.shape[1]
    _X = ((X - np.min(X, axis=0, keepdims=True)) /
          (np.max(X, axis=0, keepdims=True) - 
           np.min(X, axis=0, keepdims=True)))

    # _X = ((X - np.mean(X, axis=0, keepdims=True)) /
    #       np.var(X, axis=0, keepdims=True))
    
    return np.dot(_X, W), X


def main():
    # parse command line options
    parser = OptionParser()
    parser.add_option("--dbFinder", dest="dbFinder_file", default=None,
        help=("Results table from run_doubletFinder.R"))
    parser.add_option("--scrublet", dest="scrublet_file", default=None,
        help=("Results table from run_scrublet.py"))
    parser.add_option("--outFile", "-o", dest="out_file", default=None,
        help=("Path for output file [default: $i/combined_table.tsv]"))
    
    (options, args) = parser.parse_args()
    
    
    scrublet_dat = pd.read_csv(options.scrublet_file, header=0, sep="\t")
    dbfinder_dat = pd.read_csv(options.dbFinder_file, header=0, sep="\t")
    dbfinder_dat.cellID = [x + "-1" for x in dbfinder_dat.cellID]
    print("cell ID match: %.2f" %(np.mean(dbfinder_dat.cellID == 
                                          scrublet_dat.cellID)))
    
    _doublet_rate = np.mean(dbfinder_dat['label_fix'] == 'Doublet')
    _cmb_score, _X = average_score(dbfinder_dat['pANN_fix'], 
                                   dbfinder_dat['pANN_sweep'],
                                   scrublet_dat['score'])
    _cutoff = np.quantile(_cmb_score, 1 - _doublet_rate)
    _cmb_label = _cmb_score > _cutoff
    
    _cutoff = np.quantile(scrublet_dat['score'], 1 - _doublet_rate)
    _scb_label = scrublet_dat['score'] > _cutoff
    
    if options.out_file is None:
        dat_path = os.path.dirname(os.path.abspath(options.scrublet_file))
        out_file = dat_path + "/combined_table.tsv"
    else:
        out_file = options.out_file
        
    fid = open(out_file, "w")
    heads = ["cellID", "dbF_score_fix", "dbF_label_fix", 
             "dbF_score_sweep", "dbF_label_sweep", "scb_score", 
             "scb_label", "combined_score", "combined_label"]
    fid.writelines("\t".join(heads) + "\n")
    for i in range(len(scrublet_dat.cellID)):
        out_list = [scrublet_dat.cellID[i], 
                    "%.3f" %dbfinder_dat['pANN_fix'][i],
                    str(dbfinder_dat['label_fix'][i] == "Doublet"),
                    "%.3f" %dbfinder_dat['pANN_sweep'][i],
                    str(dbfinder_dat['label_sweep'][i] == "Doublet"),
                    "%.3f" %scrublet_dat['score'][i], str(_scb_label[i]),
                    "%.3f" %_cmb_score[i], str(_cmb_label[i])]
                    
        fid.writelines("\t".join(out_list) + "\n")
    fid.close()
    
    
if __name__ == "__main__":
    main()
    