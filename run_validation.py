
"""
Compares :
 - DIMet, against
 - pure statistical ranksum test
both used Benjamini Hochberg  correction

this script uses the output of simul_cmd.sh

"""
#
import numpy as np
from scipy.stats import contingency
from sklearn.metrics import confusion_matrix, f1_score
import numpy as np
import pandas as pd

from run_analysis import descri_args


def evaluate_method(ground_truth, predictions, method, sample_size):
    # Confusion matrix
    tn, fp, fn, tp = confusion_matrix(ground_truth, predictions).ravel()

    # Metrics
    recall = tp / (tp + fn)  # Sensitivity
    precision = tp / (tp + fp)
    specificity = tn / (tn + fp)

    f1 = f1_score(ground_truth, predictions)

    return {
        "Method": method,
        "Sample size": sample_size,
        "False positives (FP)": fp,
        "True negatives (TN)": tn,

        "Recall": np.around(recall, 2),
        "Precision": np.around(precision,2),
        "Specificity": np.around(specificity, 2),
        "F1 score": np.around(f1, 2)
    }
    # Evaluate both methods


def print_initial_check(df)->None:
    print("true positives:",
          df.loc[df['ground_truth_biomarker'] == 1, :].shape[0])
    print("true negatives:",
          df.loc[df['ground_truth_biomarker'] == 0, :].shape[0])

    print("dimet positives:", df.loc[df['dimet_biomarker'] == 1, :].shape[0])
    print("dimet negatives:", df.loc[df['dimet_biomarker'] == 0, :].shape[0])

    print("stat alone positives:",
          df.loc[df['test_alone_biomarker'] == 1, :].shape[0])
    print("stat alone negatives:",
          df.loc[df['test_alone_biomarker'] == 0, :].shape[0])


def add_response_columns(df):
    df['ground_truth_biomarker'] = 0
    df.loc[(df['d_over_s'] >= 0), 'ground_truth_biomarker'] = 1

    df['dimet_biomarker'] = 0
    df.loc[(df['d_over_s'] >= -0.2) & (
                df['padj'] <= 0.05), 'dimet_biomarker'] = 1

    df['test_alone_biomarker'] = 0
    df.loc[df['padj'] <= 0.05, 'test_alone_biomarker'] = 1

    return df

if __name__ == "__main__":

    print("Note: verify that 'group_sizes' variable has the same ")
    groups_sizes = [20, 15, 10, 7, 6, 5, 4]
    list_of_frames_quality = list()
    for gr_size in groups_sizes:
        file_path = f"../output/result_dsta-l_m1000-overlap500-a{gr_size}-b{gr_size}_10-999.tsv"
        assert "result" in file_path, "Error, 'result' must be a prefix of the input file name"
        df = pd.read_csv(file_path, sep='\t', index_col=0)
        print(df.head())

        df = add_response_columns(df)
        print_initial_check(df)

        ## compare dimet against pure statistical alone
        results = [
            evaluate_method(df['ground_truth_biomarker'].to_numpy(),
                            df['test_alone_biomarker'].to_numpy(),
                            'standard', sample_size=str(gr_size)),
            evaluate_method(df['ground_truth_biomarker'].to_numpy(),
                            df['dimet_biomarker'].to_numpy(),
                            "DIMet", sample_size=str(gr_size))
        ]

        df = pd.DataFrame(results)

        #print(df[["Method", "False sositives (FP)",
         #         "Recall", "Precision", "Specificity", "F1 score"]])

        list_of_frames_quality.append(df)

    out_df = pd.concat(list_of_frames_quality)
    out_df.to_csv("../output/validation-dimet.tsv", sep='\t')
    #file_name_ok = tabl.split("/")[-1].replace(".tsv", "")
    #output_dir = args.out_dir








