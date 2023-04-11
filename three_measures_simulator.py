"""
Data simulator with m variables, n subjects representing a,b groups
 and for each a and b groups card_a card_b respectively


Note : works with dimet conda environment

example:
python3 data_simulator.py --m 19 --m_with_overlap 10 --card_a 6 --card_b 7 \
      --value_min_abundances_df 90 --value_max_abundances_df 777 $OUTPUT_DIR

"""
import os
import random
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import argparse


def simul_data_args():
    # prog="python -m data_simulator",
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('out_dir', type=str,
                        help="directory to write results, in absolute path")

    parser.add_argument("--m", type=int, help="total number of desired variables")

    parser.add_argument("--m_with_overlap", type=int, help="number of desired variables that are overlapping\
                    (have negative distance). Must be equal or ")

    parser.add_argument( "--card_a",  type=int,
                        help="number of subjects in group a")

    parser.add_argument("--card_b", type=int,
                        help="number of subjects in group b")

    parser.add_argument("--timepoints", type=int,
                        help="how many timepoints? example 3,4,5,6...")

    parser.add_argument("--value_min_abundances_df", type=float,
                        help="minimal value for entire dataset")

    parser.add_argument("--value_max_abundances_df", type=float,
                        help="maximal value for entire dataset")

    parser.add_argument("--distances_simul_method", default="dsta-l", type=str,
                        help="[dsta-l | dsta-u | dsta-t]\n\
                         \t   Respectively, linear interpolation, uniform distribution or triangular distribution, to simulate distances.\
                         -Note: this parameter does not affect simulation of the values for each a and b groups (allways uniform)-")

    parser.add_argument('--plot_data', action=argparse.BooleanOptionalAction, default=True)

    return parser


def generate_ranges_tuples(value_min_abundances_df, value_max_abundances_df, m) -> list:
    """
    generates a list of tuples, each being a range
    example output : [(139.11, 1175.4), (85.042, 1085.22), (714.08, 1169.5), ...
    length is m
    """
    total_plage = value_max_abundances_df - value_min_abundances_df
    #critical_max_span = total_plage #* 0.9  # no more than 90% of plage
    critical_min_span = total_plage * 0.01  # no less than 1% of plage
    span_list = np.linspace(critical_min_span, total_plage, num=m)
    np.random.shuffle(span_list)
    # generate each time two numbers
    out_l = list()
    k = 0
    while len(out_l) < m:
        #random_span = np.random.uniform(critical_min_span, critical_max_span, 1)[0]
        x_inf_extreme_possible = value_max_abundances_df - span_list[k]
        random_initial_point = np.random.uniform(value_min_abundances_df, x_inf_extreme_possible, 1)[0]
        final_point = random_initial_point + span_list[k]
        # then check if they have a reasonable good span (critical max, critical min)
        out_l.append(tuple([random_initial_point, final_point]))
        k += 1
    return out_l


def distances_grouped_ranges_random_uniform(grouped_ranges_dict: dict) -> dict:
    """
    first group is Pos,
         second is Neg
    """
    odi = {'Pos': [], 'Neg': []}

    d_positive_tmp = list()
    for tup_range in grouped_ranges_dict['Pos']:
        span = abs(tup_range[0] - tup_range[1])
        d_positive = np.random.uniform(0, span, size=1)[0]
        d_positive_tmp.append(d_positive)
    # if no zeros, replace at random one value by zero (zero distance, introduce)
    if not (np.isin(0, d_positive_tmp)):
        random_index = np.random.randint(0, len(d_positive_tmp), size=1)[0]
        d_positive_tmp[random_index] = 0

    odi['Pos'] = d_positive_tmp

    d_negative_tmp = list()
    for tup_range in grouped_ranges_dict['Neg']:
        span = abs(tup_range[0] - tup_range[1])
        d_negative = np.random.uniform(-span, 0, size=1)[0] # random.uniform excludes the superior (0 here)
        d_negative_tmp.append(d_negative)
    # if no zeros, replace at random one value by zero (zero distance, introduce)

    odi['Neg'] = d_negative_tmp

    return odi


def distances_grouped_ranges_random_triang(grouped_ranges_dict: dict) -> dict:
    """
    first group is Pos,
         second is Neg
    """
    odi = {'Pos': [], 'Neg': []}

    d_positive_tmp = list()

    for tup_range in grouped_ranges_dict['Pos']:
        span = abs(tup_range[0] - tup_range[1])
        d_positive_any = np.random.triangular(left=0, mode=0, right=span, size=1)[0]
        d_positive_small = np.random.triangular(left=0, mode=0, right=span*0.1, size=1)[0]
        pick_which = np.random.randint(0,2, size=1)[0]
        if pick_which == 0:
            d_positive_tmp.append(d_positive_any)
        else:
            d_positive_tmp.append(d_positive_small)
    # if no zeros, replace at random one value by zero (zero distance, introduce)
    if not (np.isin(0, d_positive_tmp)):
        random_index = np.random.randint(0, len(d_positive_tmp), size=1)[0]
        d_positive_tmp[random_index] = 0

    odi['Pos'] = d_positive_tmp

    d_negative_tmp = list()
    for tup_range in grouped_ranges_dict['Neg']:
        span = abs(tup_range[0] - tup_range[1])
        d_negative = np.random.triangular(left=-span, right=0, mode=-span*0.25, size=1)[0] # random.uniform excludes the superior (0 here)
        d_negative_tmp.append(d_negative)
    # if no zeros, replace at random one value by zero (zero distance, introduce)

    odi['Neg'] = d_negative_tmp

    return odi


def distances_grouped_ranges_linear(grouped_ranges_dict):
    """
       first group is Pos,
            second is Neg
    """
    odi = {'Pos': [], 'Neg': []}

    d_over_s__pol = np.linspace(0, 1, num=len(grouped_ranges_dict['Pos']))
    d_over_s__nel = np.linspace(-1, -0.1, num=len(grouped_ranges_dict['Neg']))
    d_positive_tmp = list()

    for i in range(len(grouped_ranges_dict['Pos'])):
        tup_range = grouped_ranges_dict['Pos'][i]
        span = abs(tup_range[0] - tup_range[1])
        d_positive_tmp.append(span * d_over_s__pol[i])
    # if no zeros, replace at random one value by zero (zero distance, introduce)
    if not (np.isin(0, d_positive_tmp)):
        random_index = np.random.randint(0, len(d_positive_tmp), size=1)[0]
        d_positive_tmp[random_index] = 0

    odi['Pos'] = d_positive_tmp

    d_negative_tmp = list()
    for i in range(len(grouped_ranges_dict['Neg'])):
        tup_range = grouped_ranges_dict['Neg'][i]
        span = abs(tup_range[0] - tup_range[1])
        d_negative_tmp.append(span * d_over_s__nel[i])

    odi['Neg'] = d_negative_tmp

    return odi


def compute_a_b_groups(inf_lim, sup_lim, card_a, card_b, mydistance) :
    """produce a and b groups, given one single distance """
    distance = mydistance
    span = sup_lim - inf_lim
    # solving for positive and negative distances
    if distance < 0 :
        h = sup_lim - (span - (distance * -1))
        a_max = np.random.uniform(h, sup_lim, size=1)[0]
        b_min = a_max - (distance * -1)
    else:
        b_min = np.random.uniform((inf_lim + distance), sup_lim, size=1)[0]  # a number between min + distance, and max
        a_max = b_min - distance

    if card_a > 2 :
        a_values = np.random.uniform(inf_lim, a_max, size=card_a-2) # the two extremes to be added next line
        a_interval = np.append(a_values, np.array([inf_lim, a_max]) , axis=None)
    else:
        a_interval = np.array([inf_lim, a_max])

    if card_b > 2:
        b_values = np.random.uniform(b_min, sup_lim, size= card_b-2) # the two extremes to be added next line
        b_interval = np.append(b_values, np.array([b_min, sup_lim]), axis=None)
    else:
        b_interval = np.array([b_min, sup_lim])

    return a_interval, b_interval


def overlap_symmetric(x: np.array, y: np.array) -> int:
    # Credits: Claire Lescoat, Macha Nikolski, Benjamin Dartigues, Cedric Usureau, Aurélien Barré, Hayssam Soueidan
    a = [np.nanmin(x), np.nanmax(x)]
    b = [np.nanmin(y), np.nanmax(y)]

    a = np.asarray(a, dtype=np.float64)
    b = np.asarray(b, dtype=np.float64)

    overlap = np.nanmax([a[0], b[0]]) - np.nanmin([a[1], b[1]])
    return overlap


def generate_values_df(grouped_ranges_dict, grouped_distances_dict, card_a, card_b):
    outl = []
    i = 1
    for type_d in ['Neg', 'Pos']:
        distances = grouped_distances_dict[type_d]
        for j in range(len(distances)):
            distance = distances[j]
            inf_lim = min(grouped_ranges_dict[type_d][j])
            sup_lim = max(grouped_ranges_dict[type_d][j])
            a_values, b_values = compute_a_b_groups(inf_lim, sup_lim, card_a, card_b, distance)

            microdf_a = pd.DataFrame({"group": ["a" for i in range(card_a)], "value": a_values})
            microdf_b = pd.DataFrame({"group": ["b" for i in range(card_b)], "value": b_values})
            microdf = pd.concat([microdf_a, microdf_b], axis=0)
            microdf.columns = ["group", "value"]
            nums_a = [str(i + 1) for i in range(card_a)]
            nums_b = [str(i + 1) for i in range(card_b)]
            microdf["str_indiv"] = ""
            microdf.loc[microdf['group'] == "a", "str_indiv"] = nums_a
            microdf.loc[microdf['group'] == "b", "str_indiv"] = nums_b
            microdf["indiv"] = microdf['group'] + "_" + microdf['str_indiv']
            microdf = microdf.drop(columns=["str_indiv"])
            microdf["distance"] = distance  # same value all rows : distance a,b
            # calcu_dist = overlap_symmetric(a_interval, b_interval)
            # microdf["computed_d"] = calcu_dist
            microdf["var"] = f'{str(i)}'
            microdf["type_dist"] = type_d
            outl.append(microdf)
            i += 1

    df_conc = pd.concat(outl, axis=0)

    return df_conc


def plot_df_two_row(df, name) -> None:
    df["group"] = pd.Categorical(df["group"].tolist(), categories=["a", "b"] )
    sns.set_style("darkgrid")
    sns.set_palette("Dark2")

    df['distance'] = df['distance'].round(3)
    df["xlab"] = ""

    def annotate(data, **kws):
        ax = plt.gca()
        choo = data.distance.unique()[0]
        ax.text(0.01, -0.05, f"d= {choo}", transform=ax.transAxes)
        # ax.axhline(y=200, color='red')
    # goldenrod
    my_pal = {"a": "darkorange", "b": "rebeccapurple"}
    my_pal_vio = {"a": "bisque", "b": "lavender"}
    # two rows, one for "type_dist"
    g = sns.FacetGrid(df, col="var",  row="type_dist", margin_titles=True,
                      height=5, aspect=0.3, sharey=False)
    g.map_dataframe(sns.violinplot, x="xlab", y="value", data=df,
                    cut=0, # this cut in true max and min values
                    hue="group", palette=my_pal_vio, linewidth=.01, inner="quartile")
    g.map_dataframe(sns.swarmplot, x="xlab", y="value", data=df, hue="group", dodge=True,
                    palette=my_pal)
    g.map_dataframe(annotate)
    g.set_xlabels("")
    g.add_legend()
    plt.savefig(name)
    plt.close()


def plot_df_onerow(df, name) -> None:
    df["group"] = pd.Categorical(df["group"].tolist(), categories=["a", "b"] )
    sns.set_style("darkgrid")
    sns.set_palette("Dark2")

    df['distance'] = df['distance'].round(3)
    df = df.sort_values('distance', ascending=True)
    df["xlab"] = ""

    def annotate(data, **kws):
        ax = plt.gca()
        choo = data.distance.unique()[0]
        ax.text(0.01, -0.05, f"d= {choo}", transform=ax.transAxes)
        # ax.axhline(y=200, color='red')
    # goldenrod
    my_pal = {"a": "darkorange", "b": "rebeccapurple"}
    my_pal_vio = {"a": "bisque", "b": "lavender"}
    g = sns.FacetGrid(df, col="var",  margin_titles=True,
                      height=5, aspect=0.3, sharey=False)
    g.map_dataframe(sns.violinplot, x="xlab", y="value", data=df,
                    cut=0, # this cut in true max and min values
                    hue="group", palette=my_pal_vio, linewidth=.01, inner = "quartile")
    g.map_dataframe(sns.swarmplot, x="xlab", y="value", data=df, hue="group", dodge=True,
                    palette=my_pal)
    g.map_dataframe(annotate)
    g.set_xlabels("")
    g.add_legend()
    plt.savefig(name)
    plt.close()


def generate_abundance(m, card_a, card_b, value_min_abundances_df,
                       value_max_abundances_df, distances_simul_method):
    print(f"Generating abundances data:\
          \nfixed parameters: m {m},  m_with_overlap {m_with_overlap}, \
          card_a={card_a}, card_b={card_b}, \
          and for the whole dataframe min,max: [{value_min_abundances_df},{value_max_abundances_df}]")

    n = card_a + card_b  # total number of subjects

    ranges_list = generate_ranges_tuples(value_min_abundances_df, value_max_abundances_df, m)

    m_not_overlap = m - m_with_overlap
    m_list = [m_not_overlap, m_with_overlap]

    grouped_ranges_dict = {'Pos' : ranges_list[:m_list[0]],
                           'Neg' : ranges_list[- m_list[1]:]}

    if distances_simul_method == "dsta-u":
        grouped_distances_dict = distances_grouped_ranges_random_uniform(grouped_ranges_dict)
    elif distances_simul_method == "dsta-t":
        grouped_distances_dict = distances_grouped_ranges_random_triang(grouped_ranges_dict)
    elif distances_simul_method == "dsta-l" :
        grouped_distances_dict = distances_grouped_ranges_linear(grouped_ranges_dict)
    else:
        "distances_simul_method not recognized"

    df = generate_values_df(grouped_ranges_dict, grouped_distances_dict, card_a, card_b)
    print("minimal value dataframe : ",  df['value'].min())
    print("maximal : ", df['value'].max())

    pre_out = df[['indiv', 'value', 'var']]
    pre_out['var'] = "var_" + pre_out['var']
    df_abundance = pd.pivot(pre_out, index='indiv', columns='var', values='value')
    return df_abundance


def generate_fraccontri(df_abundance, metada_df):
    df_abundance




    return df


# for visualization only: isotopologues
def generate_uptolabel10_equal_isotopologues(df_abundance, metada_df):
    thecolumns = []
    return 0




if __name__ == "__main__":
    parser = simul_data_args()
    args = parser.parse_args()
    m = args.m # variables (metabolites or whatever)
    assert m <= 40, "too many variables, max 40"
    m_with_overlap = args.m_with_overlap
    card_a = args.card_a # group a size
    card_b = args.card_b# group b size
    value_min_abundances_df = args.value_min_abundances_df
    value_max_abundances_df = args.value_max_abundances_df
    distances_simul_method = args.distances_simul_method
    plot_data = args.plot_data
    out_dir = args.out_dir

    # m = 15 # variables (metabolites or whatever)
    # m_with_overlap = 5
    # card_a = 12  # group a size
    # card_b = 11 # group b size
    # value_min_abundances_df = 65  # do not put zero ?
    # value_max_abundances_df = 1110
    # distances_simul_method = "dsta-l" #  dsta-u, or dsta-l or dsta-t
    # plot_data = True

    assert m_with_overlap <= m, "Error, the m_with_overlap cannot be superior to m"
    assert value_min_abundances_df < value_max_abundances_df, "Error, the min value is not inferior to the max"
    r = distances_simul_method

    os.chdir(out_dir)
    if not os.path.exists("simulated_data/"):
        os.makedirs("simulated_data/")

    if not os.path.exists("figures/"):
        os.makedirs("figures/")

    terv = f"{value_min_abundances_df}-{value_max_abundances_df}"
    name_file_out = f"simulated_data/data_{r}_m{m}-overlap{m_with_overlap}-a{card_a}-b{card_b}_{terv}.tsv"
    plotname = f"figures/{r}_m{m}-overlap{m_with_overlap}-a{card_a}-b{card_b}_{terv}.pdf"

    df_list = list()
    for ti in range(args.timepoints):
        print("abundances, simulating timepoint : ", ti )
        df_abundance = generate_abundance(m, card_a, card_b,
                                       value_min_abundances_df,
                                        value_max_abundances_df,
                                      distances_simul_method)
        preindex = df_abundance.index

        tups = [list(i.split("_")) for i in preindex]
        indexnew = list()
        for tu in tups:
            indexnew.append(tu[0]+"_"+str(ti) + "-" + tu[1])


        df_abundance.index = indexnew
        df_list.append(df_abundance)

    df_abundance = pd.concat(df_list)
    print(df_abundance.shape)

    print("doing the metadata")
    individs = df_abundance.index
    metada_df = pd.DataFrame({'former_name': individs,
                              'sample': individs,
                              'condition': [i[0] for i in individs],
                              'timepoint': [i[2] for i in individs],
                              'timenum': [i[2] for i in individs],
                              'short_comp': ['en' for i in individs]
                              })
    print(metada_df)
    file_metada = name_file_out.split("/")[-1].\
        replace(".tsv", ".csv").replace("data_", "meta_")
    metada_df.to_csv(file_metada, index=False)
    df_abundance.to_csv(name_file_out, sep='\t', index=True)
    print("\n  Done Abundance")

    print('\nGenerating synthetic fractional contributions')

    df_fraccontri = generate_fraccontri(df_abundance,
                                        metada_df)

    # for visualization only: isotopologues
    df_isotopol = generate_uptolabel10_equal_isotopologues(m, card_a, card_b,
                                                  metada_df)















