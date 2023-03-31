"""
plot different metrics from simulateddata
"""
import os
import numpy as np
import pandas as pd
import argparse


def metada_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('in_file', type=str,
                        help="data file to describe, absolute path")

    parser.add_argument('dir_out', type=str,
                        help="directory to write results, in absolute path")

    return parser

if __name__ == "__main__":

    parser = metada_args()
    args = parser.parse_args()

    tabl = args.in_file
    df = pd.read_csv(tabl, sep='\t', index_col=0)

    individs = df.index

    metada_df = pd.DataFrame({'former_name': individs,
                              'sample': individs,
                              'condition': [i[0] for i in individs],
                              'timepoint': ['0' for i in range(len(individs))],
                              'timenum': [0 for i in range(len(individs))]})

    file_out = args.dir_out + tabl.split("/")[-1].replace(".tsv", ".csv").replace("data_", "meta_")
    metada_df.to_csv(file_out, index=False)