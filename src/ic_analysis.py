import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from pathlib import Path
import os
from math import ceil

def load_ic_data(input_folder):
    df_dict = {}
    ic_path = Path(input_folder)

    # Iterowanie przez pliki CSV i tworzenie DataFrame'ów
    for csv_file in ic_path.glob("*.csv"):
        network, dataset = csv_file.stem.split("_")[:2]
        if network not in df_dict:
            df_dict[network] = {}
        df_dict[network][dataset] = pd.read_csv(csv_file)

    return df_dict

def prepare_dataframe(df_dict):
    df_list = []

    # Przetwarzanie i agregowanie danych
    for network in df_dict:
        for dataset, dataset_df in df_dict[network].items():
            dataset_df["species_network"] = network
            dataset_df["species"] = network[:4]  # Zakładam, że pierwsze 4 litery to gatunek
            dataset_df["network"] = network[5:]  # Zakładam, że reszta po _ to sieć
            dataset_df["dataset"] = dataset
            df_list.append(dataset_df)

    df_all = pd.concat(df_list, ignore_index=True)
    return df_all

def plot_ic_distributions(df_all, output_folder):
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    ic_max = ceil(df_all['IC'].max() * 2) / 2

    g = sns.FacetGrid(df_all, row="species", col="network", height=4, aspect=1.5)

    def overlay_cdf(data, color='k', **kwargs):
        ax2 = plt.gca().twinx()
        sns.ecdfplot(data=data, ax=ax2, color=color, alpha=0.9)
        ax2.set_ylabel('CDF', color=color)
        ax2.tick_params(axis='y', labelcolor=color)
        ax2.grid(True)

    g.map(sns.histplot, 'IC', binwidth=0.5, binrange=(0, ic_max))
    g.map(overlay_cdf, 'IC', color='red')

    g.fig.subplots_adjust(top=0.94)
    g.fig.suptitle("Histograms and CDF plots for 6 networks | 5mers")
    
    hist_cdf_plot_path = os.path.join(output_folder, "histograms_cdf.png")
    g.savefig(hist_cdf_plot_path)
    plt.close(g.fig)

def plot_ic_density(df_all, output_folder):
    sns.displot(df_all, x="IC", hue="species_network", kind="kde")
    plt.grid()
    plt.title("Density plots overlayed for 6 networks | 5mers")
    
    density_plot_path = os.path.join(output_folder, "density_plots.png")
    plt.savefig(density_plot_path)
    plt.close()

def plot_single_histograms(df_all, output_folder):
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    
    for species_network, data in df_all.groupby('species_network'):
        plt.figure(figsize=(8, 6))
        sns.histplot(data['IC'], bins=30, kde=False)
        plt.title(f"IC Histogram for {species_network}")
        plt.xlabel('IC')
        plt.ylabel('Count')
        plt.grid()

        single_histogram_path = os.path.join(output_folder, f"{species_network}_histogram.png")
        plt.savefig(single_histogram_path)
        plt.close()

def ic_plots(input_folder, output_folder):
    df_dict = load_ic_data(input_folder)
    df_all = prepare_dataframe(df_dict)
    
    plot_ic_distributions(df_all, output_folder)
    plot_ic_density(df_all, output_folder)
    plot_single_histograms(df_all, output_folder)

