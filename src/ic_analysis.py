import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from pathlib import Path
import os
from math import ceil

def load_ic_data(input_folder):
    """
    Load Information Content (IC) data from CSV files in the specified folder.

    Args:
        input_folder (str): The folder where the CSV files are located.

    Returns:
        pd.DataFrame: A combined DataFrame with IC data and additional metadata columns.
    """
    df_list = []
    ic_path = Path(input_folder)

    # Iterating through CSV files and creating DataFrames
    for csv_file in ic_path.glob("*.csv"):
        filename_parts = csv_file.stem.split("_")
        if len(filename_parts) >= 3:
            dataset, network, class_name = filename_parts[:3]
        else:
            print(f"Filename {csv_file.name} does not match expected pattern.")
            continue

        df = pd.read_csv(csv_file)
        df["dataset"] = dataset
        df["network"] = network
        df["class"] = class_name
        df_list.append(df)

    if df_list:
        df_all = pd.concat(df_list, ignore_index=True)
    else:
        df_all = pd.DataFrame()

    return df_all

def prepare_dataframe(df_all):
    """
    Prepare the DataFrame by adding any additional columns or processing.

    Args:
        df_all (pd.DataFrame): The combined DataFrame with IC data.

    Returns:
        pd.DataFrame: The processed DataFrame.
    """
    # No additional processing needed based on current requirements
    return df_all

def plot_ic_distributions(df_all, output_folder, k):
    """
    Generate histograms and CDF plots for IC values across datasets and networks.

    Args:
        df_all (pd.DataFrame): The combined DataFrame with IC data.
        output_folder (str): The folder where the output plots will be saved.
        k (int): The k-mer length used in IC calculation.
    """
    import matplotlib.pyplot as plt
    import seaborn as sns
    from math import ceil
    import os

    if df_all.empty:
        print("No data to plot.")
        return

    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    ic_max = ceil(df_all['IC'].max() * 2) / 2  # Ceiling to the nearest 0.5

    g = sns.FacetGrid(df_all, row="dataset", col="network", height=4, aspect=1.5)

    # Overlay CDF plot on top of the histograms
    def overlay_cdf(data, **kwargs):
        ax = plt.gca()
        ax2 = ax.twinx()
        sns.ecdfplot(data=data, ax=ax2, color='red', alpha=0.9)
        ax2.set_ylabel('CDF', color='red')
        ax2.tick_params(axis='y', labelcolor='red')
        ax2.grid(True)

    g.map(sns.histplot, 'IC', binwidth=0.5, binrange=(0, ic_max), color='blue')
    g.map(overlay_cdf, 'IC')

    # Adjust titles to prevent overlapping
    g.set_titles(col_template="{col_name}", row_template="{row_name}", size=12)
    g.fig.subplots_adjust(top=0.9)  # Increase top margin to make space for suptitle
    g.fig.suptitle(f"Histograms and CDF Plots | {k}-mers", fontsize=16, y=0.98)

    hist_cdf_plot_path = os.path.join(output_folder, "histograms_cdf.png")
    plt.savefig(hist_cdf_plot_path)
    plt.close(g.fig)


def plot_ic_density(df_all, output_folder, k):
    """
    Generate kernel density estimates (KDE) plots for IC values across classes.

    Args:
        df_all (pd.DataFrame): The combined DataFrame with IC data.
        output_folder (str): The folder where the output plots will be saved.
        k (int): The k-mer length used in IC calculation.
    """
    if df_all.empty:
        print("No data to plot.")
        return

    g = sns.displot(df_all, x="IC", hue="class", kind="kde", height=5, aspect=1.5)
    # Update the title
    g.fig.suptitle(f"Density plots overlaid | {k}-mers", y=0.95)
    # Adjust top to prevent the title from being cut off
    g.fig.subplots_adjust(top=0.9)
    # Add grid to each subplot
    for ax in g.axes.flatten():
        ax.grid(True)
    
    density_plot_path = os.path.join(output_folder, "density_plots.png")
    g.savefig(density_plot_path)
    plt.close()

def plot_single_histograms(df_all, output_folder):
    """
    Generate individual histograms for IC values per dataset_network_class.

    Args:
        df_all (pd.DataFrame): The combined DataFrame with IC data.
        output_folder (str): The folder where the output histograms will be saved.
    """
    if df_all.empty:
        print("No data to plot.")
        return

    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    
    df_all['dataset_network_class'] = df_all['dataset'] + '_' + df_all['network'] + '_' + df_all['class']

    # Generate a histogram for each dataset_network_class
    for group_name, data in df_all.groupby('dataset_network_class'):
        plt.figure(figsize=(8, 6))
        sns.histplot(data['IC'], bins=30, kde=False)
        plt.title(f"IC Histogram for {group_name}")
        plt.xlabel('IC')
        plt.ylabel('Count')
        plt.grid()

        # Sanitize group_name for filename
        sanitized_group_name = group_name.replace('/', '_')
        single_histogram_path = os.path.join(output_folder, f"{sanitized_group_name}_histogram.png")
        plt.savefig(single_histogram_path)
        plt.close()

def ic_plots(input_folder, output_folder, k):
    """
    Main function to generate all IC plots: histograms, density plots, and single histograms.

    Args:
        input_folder (str): The folder containing the IC CSV files.
        output_folder (str): The folder where the generated plots will be saved.
        k (int): The k-mer length used in IC calculation.
    """
    df_all = load_ic_data(input_folder)
    if df_all.empty:
        print("No data loaded. Exiting.")
        return

    df_all = prepare_dataframe(df_all)
    
    plot_ic_distributions(df_all, output_folder, k)
    plot_ic_density(df_all, output_folder, k)
    plot_single_histograms(df_all, output_folder)

def compute_default_ic_threshold(k):
    """Compute the default IC threshold as (k + 1) / 2."""
    return (k + 1) / 2
