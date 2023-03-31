"""
plot different metrics from simulateddata
"""
import os
import argparse
#import random
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
import locale
from scipy import stats

from data_simulator import plot_df_onerow, overlap_symmetric


def descri_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('in_file', type=str,
                        help="data file to describe, absolute path")

    parser.add_argument('out_dir', type=str,
                        help="directory to write results, in absolute path")

    parser.add_argument('--redo_plot_orig_data', action=argparse.BooleanOptionalAction, default=False)
    return parser


def compute_gmean_nonan(anarray):
    anarray = np.array(anarray, dtype=float)
    anarray = anarray[~np.isnan(anarray)]
    if sum(anarray) == 0:  # replicates all zero
        outval = 0
    else:
        outval = stats.gmean(anarray)
    return outval


def compute_intervals_relation(x ,y):
    local_span_x = max(x) - min(x)
    local_span_y = max(y) - min(y)
    lsx = local_span_x
    lsy = local_span_y
    times_the_bigger_in_the_smaller = ( min(lsx ,lsy) +1) / (max(lsx, lsy) + 1)
    return times_the_bigger_in_the_smaller


def compute_distance_over_span(d, combined_a_b_values):
    span = max(combined_a_b_values) - min(combined_a_b_values)
    return float(d / span)


def compute_b_versus_a_fold_change(b_values: np.array, a_values: np.array):
    fc = compute_gmean_nonan(b_values) / compute_gmean_nonan(a_values)
    return fc


def compute_abs_normalized_diff(b_values: np.array, a_values: np.array):
    m_a = compute_gmean_nonan(b_values)
    m_b = compute_gmean_nonan(a_values)
    denom = m_a + m_b
    result = (m_a - m_b) / denom
    return abs(result)


def compute_ranksums_allH0(groupB: np.array, groupA: np.array):
    # The Wilcoxon rank-sum test tests the null hypothesis that two sets of measurements are drawn from the same distribution
    # ‘two-sided’: one of the distributions (underlying x or y) is stochastically greater than the other.
    # ‘less’: the distribution underlying x is stochastically less than the distribution underlying y.
    #  ‘greater’: the distribution underlying x is stochastically greater than the distribution underlying y.
    groupA = groupA[~np.isnan(groupA)]
    groupB = groupB[~np.isnan(groupB)]
    if np.median(groupB) > np.median(groupA):
        detected_alternative = 'greater'
    elif np.median(groupB) < np.median(groupA):
        detected_alternative = 'less'
    else:
        detected_alternative = 'two-sided'
    stat, pval = stats.ranksums(groupB, groupA, alternative=detected_alternative)
    return stat, pval


def do_descriptor_df(df, groups):
    m = df['var_num'].max()
    intervals_rel = []
    distance = []
    d_over_s = []
    fold_change_b_a = []
    absolute_diff = []
    pvalue_pkg = []
    for i in range(1, (m + 1)):
        sub_df = df.loc[df['var_num'] == i, :]
        groups_both_vals = sub_df['value'].to_numpy()
        groupA = sub_df.loc[sub_df['group'] == groups[0], 'value'].to_numpy()
        groupB = sub_df.loc[sub_df['group'] == groups[1], 'value'].to_numpy()
        len_grA = len(groupA)
        len_grB = len(groupB)
        # comparison is B/ A
        distance_here = overlap_symmetric(groupA, groupB)
        distance.append(distance_here)
        intervals_rel.append(compute_intervals_relation(groupB, groupA))
        d_over_s.append(compute_distance_over_span(distance_here, groups_both_vals))
        fold_change_b_a.append(compute_b_versus_a_fold_change(groupB, groupA))
        absolute_diff.append(compute_abs_normalized_diff(groupB, groupA))
        stat, pval = compute_ranksums_allH0(groupB, groupA)
        pvalue_pkg.append(pval)

    df_out_pkg = pd.DataFrame({"intervals_rel": intervals_rel,
                               "distance": distance,
                               "d_over_s": d_over_s,
                               "fold_change_b_a": fold_change_b_a,
                               "absolute_diff": absolute_diff,
                               "pvalue_pkg": pvalue_pkg})

    df_out_pkg.index = ["var_" + str(i) for i in range(1, (m + 1))]
    df_out_pkg['var'] = df_out_pkg.index

    return df_out_pkg, len_grA, len_grB



def give_reduced_df( df, ddof ):
    rownames = df.index
    df.index = range(len(rownames))  # index must be numeric because compute reduction accepted
    df_red = compute_reduction(df, ddof)  # reduce
    df_red.index = rownames
    return df_red


def compute_reduction(df, ddof):
    """
    modified, original from ProteomiX
    johaGL 2022: if all row is zeroes, set same protein_values
    """
    res = df.copy()
    for protein in df.index.values:
        # get array with abundances values
        protein_values = np.array(
            df.iloc[protein].map(lambda x: locale.atof(x) if type(x) == str else x) )
        # return array with each value divided by standard deviation, row-wise
        if sum(protein_values) == 0:
            reduced_abundances = protein_values  # because all row is zeroes
        else:
            reduced_abundances = protein_values / np.nanstd(protein_values, ddof=ddof)

        # replace values in result df
        res.loc[protein] = reduced_abundances
    return res

if __name__ == "__main__":

    parser = descri_args()
    args = parser.parse_args()
    tabl = args.in_file
        #"../simulated_data/data_dsta-l_m50-overlap45-a11-b12_67.0-1335.0.tsv"
    #tabl =  "../simulated_data/data_unif_m1000-overlap300-a33-b32_0-1.tsv"
    dir_out = args.dir_out + tabl.split("/")[-1].replace(".tsv", "").replace("data_", "") + "/"
    if not os.path.exists(dir_out):
        os.makedirs(dir_out)

    redo_plot_orig_data = args.redo_plot_orig_data

    df = pd.read_csv(tabl,  sep='\t', index_col=0)
    dfredu_transposed = give_reduced_df(df.T, 0)
    dfredu = dfredu_transposed.T
    df = dfredu.copy() # ok, reduced
    df["indiv"] = df.index
    df = pd.melt(df, id_vars="indiv", var_name="var", value_name='value')
    df = df.assign(group=df['indiv'].str[0]) # the first character of indiv: a_10 ==> a

    # in all data, replace zero by the minimal of the data
    min_all_data = df['value'].min()
    df.loc[df['value'] == 0, 'value'] = min_all_data

    df = df.assign(var_num=df['var'].str.replace("var_", "").astype(int))

    m = df['var_num'].max()
    print("nb of variables: ", m)
    df = df.sort_values("var_num", ascending=True)
    # remind, var is the variable (metabolite, etc)
    # by var and by group, obtain descriptive metrics:
    groups = ['a', 'b']
    df_out_pkg, len_grA, len_grB = do_descriptor_df(df, groups)


    if redo_plot_orig_data:
        if m > 60:
            print("Won t print plot, too many variables: m= (m)")
        else : # no more than 60 variables for this kind of plot
            print("\nPrinting plot of simulated variables...")
            df_out_pkg['var'] = df_out_pkg.index
            df_sele_cols = df_out_pkg[['var', 'distance']]
            df_pl = pd.merge(df_sele_cols, df, on="var", how="left")
            sns.set_style("darkgrid")
            sns.set_palette("Dark2")
            df_pl = df_pl.sort_values("distance", ascending=True)
            plot_df_onerow(df_pl, dir_out + "1.pdf")


    #sns.set_theme(style="white")
    sns.set_theme(style="ticks")
    sns.jointplot(data=df_out_pkg,
                  x='absolute_diff', y= 'd_over_s', color ="#4CB391")
    plt.savefig(dir_out + "2.pdf")
    #plt.show()

    sns.jointplot(data=df_out_pkg,
                  x='intervals_rel', y='d_over_s', color="#4CB391")
    plt.savefig(dir_out + "3.pdf")

    #plt.show()

    sns.jointplot(data=df_out_pkg,
                  x='absolute_diff', y='intervals_rel', color="#4CB391")
    plt.savefig(dir_out + "4.pdf")
    #plt.show()

    print()

    ####################################
    ### bubbles plot pvalue viridis
    ####################################

    sns.set_theme(style="darkgrid")

    plt.figure(figsize=(10, 8))
    plt.scatter(data=df_out_pkg, x='absolute_diff', y='d_over_s', alpha=0.75,
                c='pvalue_pkg', s=(df_out_pkg['intervals_rel'].to_numpy()*100) ** 1.25, cmap="viridis_r")
    #s = (intervals_rel ** 2) * 60,
    ax = plt.gca()
    min_bub = df_out_pkg['intervals_rel'].min()
    max_bub = df_out_pkg['intervals_rel'].max()

    plt.colorbar(label="pval_non_parametric_test")
    plt.xlabel("absolute_diff")
    plt.ylabel("d_over_s")

    # make a legend:
    pws = np.linspace(min_bub, max_bub, 4).round(2)
    for pw in pws:
       #plt.scatter([], [], s=(pw ** 2) * 60, c="k", label=str(pw))
        plt.scatter([], [], s=(pw*100)**1.25, c="darkgray", label=str(pw))
    h, l = plt.gca().get_legend_handles_labels()
    plt.legend(h[1:], l[1:], labelspacing=1.2, title="intervals_rel", borderpad=1,
               frameon=True, framealpha=0.6, edgecolor="k", facecolor="w")

    plt.title(f"Comparing group B (n={len_grB}) vs A (n={len_grA}) across {m} variables. ")

    #plt.show()
    plt.savefig(dir_out + "5.pdf")
    plt.close()
    ###

    ###########
    ### same plot but classify into significant and not significant
    ###########

    df_out_pkg['p_class'] = ''
    df_out_pkg.loc[df_out_pkg['pvalue_pkg'] <= 0.05, 'p_class'] = '<= 0.05 (significant)'
    df_out_pkg.loc[df_out_pkg['pvalue_pkg'] > 0.05, 'p_class'] = '> 0.05'

    mypal_categ_pval = {'<= 0.05 (significant)': 'coral',
                        '> 0.05': 'dodgerblue'}

    sns.set_theme(style="darkgrid")

    plt.figure(figsize=(10, 8))
    plt.scatter(data=df_out_pkg, x='absolute_diff', y='d_over_s', alpha=0.5,
                c=df_out_pkg['p_class'].map(mypal_categ_pval),
                s=(df_out_pkg['intervals_rel'].to_numpy()*100) ** 1.25)

    ax = plt.gca()
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])

    # make a legend:
    min_bub = df_out_pkg['intervals_rel'].min()
    max_bub = df_out_pkg['intervals_rel'].max()
    pws = np.linspace(min_bub, max_bub, 4).round(2)
    for pw in pws:
       #plt.scatter([], [], s=(pw ** 2) * 60, c="k", label=str(pw))
        plt.scatter([], [], s=(pw*100)**1.25, c="darkgray", label=str(pw))
    h, l = plt.gca().get_legend_handles_labels()
    plt.legend(h[1:], l[1:], labelspacing=1.2, title="intervals_rel", borderpad=1,
               frameon=True, framealpha=0.6,   edgecolor="k", facecolor="w")
    patch1 = mpatches.Patch(color='coral', label='<= 0.05 (significant)')
    patch2 = mpatches.Patch(color='dodgerblue', label='> 0.05')
    plt.legend(handles=[patch1, patch2], bbox_to_anchor=(1, 0.5), loc='upper left')
    # plt.tight_layout()
    plt.title(f"Comparing group B (n={len_grB}) vs A (n={len_grA}) across {m} variables. ")

    # add dots annotations:
    stytext = dict(size=8, color='gray')
    for i, r in df_out_pkg.iterrows():
        ax.text(x=df_out_pkg.loc[i, 'absolute_diff'],
                y=df_out_pkg.loc[i, 'd_over_s'],
                s=df_out_pkg.loc[i, 'var'], **stytext)
    #plt.show()
    plt.xlabel("absolute_diff")
    plt.ylabel("d_over_s")
    plt.savefig(dir_out + "6.pdf")
    plt.close()



    def relplot_old():
        df_out_pkg['p_class'] = ''
        df_out_pkg.loc[df_out_pkg['pvalue_pkg'] <= 0.05, 'p_class'] = '<= 0.05 (significant)'
        df_out_pkg.loc[df_out_pkg['pvalue_pkg'] > 0.05, 'p_class'] = '> 0.05'

        mypal_categ_pval = {'<= 0.05 (significant)' : 'coral',
                                '> 0.05' : 'dodgerblue' }

        plt.figure(figsize=(13, 9))

        relp = sns.relplot(x="absolute_diff", y="d_over_s", hue="p_class", size='intervals_rel',
                    sizes=(25, 250),
                    alpha=.5,  palette=mypal_categ_pval, #palette="muted",
                    #height=6,
                           data=df_out_pkg)
        plt.title(f"Comparing group B (n={len_grB}) vs A (n={len_grA}) across {m} variables. ")
        ax = relp.axes[0,0]
        # add dots annotations:
        stytext = dict(size=8, color='gray')
        for i,r in df_out_pkg.iterrows():
            ax.text(x=df_out_pkg.loc[i, 'absolute_diff'],
                     y=df_out_pkg.loc[i, 'd_over_s'],
                     s=df_out_pkg.loc[i, 'var'], **stytext)

        # plt.show()
        plt.savefig(dir_out + "6-old.pdf")


    ###
    print("cases where distance is negative but padj is signif : ")
    whois_lodiff_lods_psignif  = df_out_pkg.loc[df_out_pkg['pvalue_pkg'] <= 0.05, :]
    whois_lodiff_lods_psignif = whois_lodiff_lods_psignif.loc[
        whois_lodiff_lods_psignif['d_over_s'] <= -0.9, : ]
    print(whois_lodiff_lods_psignif)

