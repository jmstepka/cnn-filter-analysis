import pandas as pd
import numpy as np
import os
import seaborn as sns
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.stats import kstest
import logomaker as lm

def load_all_network_df(directory="output_kl"):
    """
    Load all network data from CSV files and concatenate into a single DataFrame.

    Args:
        directory (str): Directory containing the KL distance CSV files.

    Returns:
        pd.DataFrame: DataFrame containing KL distances for all networks.
    """
    dfs = []
    for file in Path(directory).glob("*.csv"):
        dfs.append(pd.read_csv(file, sep="\t"))
    df = pd.concat(dfs, ignore_index=True)
    df["KL distance"] = -1 * df["KL distance"]
    return df

def tag_family(df, path, name):
    """
    Tag motifs belonging to a specific family based on a file of family TFs.

    Args:
        df (pd.DataFrame): DataFrame containing KL distance data.
        path (str): Path to the file containing family TFs.
        name (str): Name of the family to tag.

    Returns:
        pd.DataFrame: DataFrame with an additional tag column indicating family membership.
    """
    family_tfs = pd.read_csv(f"{path}", sep="\t")
    family_tfs[f"{name}_tag"] = 1
    family_df = family_tfs[["Model", f"{name}_tag"]]
    tagged = df.merge(family_df, left_on="Target", right_on="Model", how="left")
    tagged[f"{name}_tag"].fillna(0, inplace=True)
    return tagged

def get_enrichment_df(sorted_df, tag_name, shortening_thresh=1.5):
    """
    Calculate enrichment data for motifs in a family at different KL distance thresholds.

    Args:
        sorted_df (pd.DataFrame): Sorted DataFrame by KL distance.
        tag_name (str): Column name indicating family membership.
        shortening_thresh (float): KL distance threshold for reducing data density.

    Returns:
        pd.DataFrame: Enrichment data for the family across KL distances.
    """
    dfs = []
    for network in sorted_df["Network"].unique():
        enrichment_data = []
        dist_data = []
        total = 0
        family_count = 0
        for r in sorted_df[sorted_df["Network"] == network].iterrows():
            r = r[1]
            family_member = r[tag_name]
            total += 1
            family_count += family_member
            enrichment_data.append(family_count / total)
            dist_data.append(r["KL distance"])
        dict_data = {"enrichment": enrichment_data, "KL distance": dist_data}
        df = pd.DataFrame(dict_data)
        df["Network"] = network
        dfs.append(df)
    enrichment_df = pd.concat(dfs, axis=0)
    enrichment_df_shortened = pd.concat(
        [enrichment_df[enrichment_df["KL distance"] <= shortening_thresh],
         enrichment_df[enrichment_df["KL distance"] > shortening_thresh][::20]])
    return enrichment_df_shortened

def plot_kl_distributions(df_all, output_dir):
    """
    Plot histograms and CDF plots for KL distances across networks and datasets.

    Args:
        df_all (pd.DataFrame): DataFrame containing KL distance data.
        output_dir (str): Directory to save the output plot.
    """
    g = sns.FacetGrid(df_all, row="Dataset", col="Network", height=4, aspect=1.5)

    def overlay_cdf(data, color='k', **kwargs):
        ax2 = plt.gca().twinx()
        sns.ecdfplot(data, color=color, alpha=0.9, ax=ax2)
        ax2.set_ylabel('CDF', color=color)
        ax2.tick_params(axis='y', labelcolor=color)
        ax2.grid(True)

    g.map(sns.histplot, 'KL distance')
    g.map(overlay_cdf, 'KL distance', color='red')
    thresh = 5.0
    g.fig.subplots_adjust(top=0.94)
    g.fig.suptitle(f"Histograms and CDF plots for all networks | KL distance threshold = {thresh}")
    
    # Save plot
    output_path = os.path.join(output_dir, "KL_distance_distribution.png")
    g.savefig(output_path)
    plt.close()

def plot_kl_distributions_with_motif_family(df_all, family, output_dir):
    """
    Plot KL distance distributions for motifs by family membership.

    Args:
        df_all (pd.DataFrame): DataFrame with KL distance and family tags.
        family (str): Family column name to differentiate in the plot.
        output_dir (str): Directory to save the output plot.
    """
    g = sns.FacetGrid(df_all, row="Dataset", col="Network", hue=family, height=4, aspect=1.5)
    g.map(sns.histplot, 'KL distance')
    thresh = 5.0
    g.fig.subplots_adjust(top=0.94)
    g.fig.suptitle(f"Histograms by family / non-family | KL distance threshold = {thresh}")
    g.add_legend()

    # Save plot
    output_path = os.path.join(output_dir, f"KL_distribution_family_{family}.png")
    g.savefig(output_path)
    plt.close()

def plot_enrichment(enrich_df, family, output_dir):
    """
    Plot family enrichment at different KL distance thresholds.

    Args:
        enrich_df (pd.DataFrame): Enrichment data for the family.
        family (str): Name of the family.
        output_dir (str): Directory to save the output plot.
    """
    sns.lineplot(data=enrich_df, x="KL distance", y="enrichment", hue="Network")
    plt.ylabel(f"{family} proportion")
    plt.title(f"{family} proportion in different KL distance thresholds")

    # Save plot
    output_path = os.path.join(output_dir, f"KL_distance_enrichment_{family}.png")
    plt.savefig(output_path)
    plt.close()

def plot_csf_with_kstext(df_all, family, output_dir):
    """
    Plot CDF and perform the Kolmogorov-Smirnov test for family vs non-family motifs.

    Args:
        df_all (pd.DataFrame): DataFrame with KL distance and family tags.
        family (str): Family column name for testing.
        output_dir (str): Directory to save the output plot.
    """
    g = sns.FacetGrid(df_all, row="Dataset", col="Network", hue=family, height=4, aspect=1.5)

    def overlay_cdf(data, color='k', **kwargs):
        ax2 = plt.gca()
        sns.ecdfplot(data, color=color, alpha=0.9, ax=ax2)
        ax2.set_ylabel('CDF')
        ax2.grid(True)

    g.map(overlay_cdf, 'KL distance')

    titles = []
    for network in df_all["Network"].unique():
        data = df_all[df_all["Network"] == network]
        pvalue = kstest(data['KL distance'], data[data[family] == 1]['KL distance'], alternative='less').pvalue
        titles.append(f"{network}, p-value: {pvalue:.3g}")

    for ax, title in zip(g.axes.flatten(), titles):
        ax.set_title(title)

    g.fig.subplots_adjust(top=0.94)
    g.fig.suptitle(f"Kolmogorov-Smirnov Test Results for Family vs Non-family | KL distance")
    g.add_legend()

    # Save plot
    output_path = os.path.join(output_dir, f"KL_distribution_KS_test_{family}.png")
    g.savefig(output_path)
    plt.close()

def perform_kl_analysis(input_dir, family_files, family_names, output_dir):
    """
    Perform KL analysis, including plotting and family enrichment analysis.

    Args:
        input_dir (str): Directory containing KL distance data.
        family_files (list): List of paths to family motif files.
        family_names (list): List of family names.
        output_dir (str): Directory to save all output plots.
    """
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    df_all = load_all_network_df(directory=input_dir)

    for family_file, family_name in zip(family_files, family_names):
        print(f"Processing family: {family_name}")

        df_all_tagged = tag_family(df_all, family_file, family_name)

        plot_kl_distributions(df_all_tagged, output_dir)

        plot_kl_distributions_with_motif_family(df_all_tagged, family=f"{family_name}_tag", output_dir=output_dir)

        sorted_df = df_all_tagged.sort_values(by=["KL distance"])
        enrichment_df = get_enrichment_df(sorted_df, f"{family_name}_tag")
        plot_enrichment(enrichment_df, family=family_name, output_dir=output_dir)

        plot_csf_with_kstext(df_all_tagged, family=f"{family_name}_tag", output_dir=output_dir)
