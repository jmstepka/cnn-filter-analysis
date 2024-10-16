import seaborn as sns
import pandas as pd
import logomaker as lm
import matplotlib.pyplot as plt

from scipy.stats import kstest

from data_processing import get_matrix, get_tf_matrix

from math import ceil

def plot_ic_distributions(df_all):
    ic_max = ceil(df_all['IC'].max() * 2) / 2

    g = sns.FacetGrid(df_all, row="species", col="network", height=4, aspect=1.5)

    def overlay_cdf(data, color='k', **kwargs):
        
        ax2 = plt.gca().twinx()

        ax2.ecdf(data, color=color, alpha=0.9)
        ax2.set_ylabel('CDF', color=color)
        ax2.tick_params(axis='y', labelcolor=color)

        ax2.grid(True)
        #ax2.yaxis.grid(True)


    g.map(sns.histplot, 'IC', binwidth=0.5, binrange=(0, ic_max))
    g.map(overlay_cdf, 'IC', color='red')

    g.fig.subplots_adjust(top=0.94)
    g.fig.suptitle("Histograms and CDF plots for 6 networks | 5mers")

def plot_ic_density(df_all):
    sns.displot(df_all, x="IC", hue="species_network", kind="kde")
    plt.grid()

    plt.title("Density plots overlayed for 6 networks | 5mers")

def plot_kl_distributions(df_all):
    g = sns.FacetGrid(df_all, row="Species", col="network_type", height=4, aspect=1.5)

    def overlay_cdf(data, color='k', **kwargs):
        
        ax2 = plt.gca().twinx()

        ax2.ecdf(data, color=color, alpha=0.9)
        ax2.set_ylabel('CDF', color=color)
        ax2.tick_params(axis='y', labelcolor=color)

        ax2.grid(True)


    g.map(sns.histplot, 'KL distance')
    g.map(overlay_cdf, 'KL distance', color='red')
    tresh = 5.0
    g.fig.subplots_adjust(top=0.94)
    g.fig.suptitle(f"Histograms and CDF plots for 6 networks | 5mers | thresh = {tresh}")

def plot_kl_distributions_with_motif_family(df_all, family):
    g = sns.FacetGrid(df_all, row="Species", col="network_type", hue=family, height=4, aspect=1.5)

    g.map(sns.histplot, 'KL distance')

    thresh = 5.0
    g.fig.subplots_adjust(top=0.94)
    g.fig.suptitle(f"Histograms by stripe / non-stripe | 5mers | thresh = {thresh}")

def logo_row(row):
    fig, axs = plt.subplots(nrows=2)
    kl = row["KL distance"]
    infomat = get_matrix(row)
    tf_infomat = get_tf_matrix(row)
    lm.Logo(infomat, ax=axs[0])
    lm.Logo(tf_infomat, ax=axs[1])
    fig.suptitle(kl,fontsize=21)

def plot_enrichment(enrich_df, family):
    sns.lineplot(enrich_df, x="KL distance", y="enrichment", hue="network")
    plt.ylabel(f"{family} proportion")
    plt.title(f"{family} proportion in different distance thresholds")

def plot_csf_with_kstext(all_df_sorted, family):
    g = sns.FacetGrid(all_df_sorted, row="Species", col="network_type", hue=family, height=4, aspect=1.5)

    def overlay_cdf(data, color='k', **kwargs):
        
        ax2 = plt.gca()

        ax2.ecdf(data, color=color, alpha=0.9)
        ax2.set_ylabel('CDF')
        ax2.tick_params(axis='y')

        ax2.grid(True)
        #ax2.yaxis.grid(True)


    #g.map(sns.histplot, 'KL distance')
    g.map(overlay_cdf, 'KL distance')


    titles = []
    for network in get_network_dataset_dict():
        data = all_df_sorted[all_df_sorted["Network"] == network]
        pvalue = kstest(data['KL distance'], data[data[family] == 1]['KL distance'], alternative='less').pvalue
        titles.append(f"{network}, {pvalue:.3g}")

    for ax, title in zip(g.axes.flatten(),titles):
        ax.set_title(title)

    g.fig.subplots_adjust(top=0.94)
    g.fig.suptitle(f"Histograms by stripe / non-stripe | 5mers | thresh = {tresh}")
    g.add_legend()