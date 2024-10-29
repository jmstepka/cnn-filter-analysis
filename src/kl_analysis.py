import pandas as pd
import numpy as np
import os
import seaborn as sns
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.stats import kstest
from statsmodels.stats.multitest import multipletests

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

def load_tf_family_data(tf_family_file):
    """
    Load the TF family data from a CSV file.

    Args:
        tf_family_file (str): Path to the TF family CSV file.

    Returns:
        pd.DataFrame: DataFrame containing TF family data.
    """
    tf_family_df = pd.read_csv(tf_family_file)
    return tf_family_df

def merge_tf_family_data(df_kl, tf_family_df):
    """
    Merge the TF family data into the KL distance dataframe.

    Args:
        df_kl (pd.DataFrame): KL distance dataframe.
        tf_family_df (pd.DataFrame): TF family dataframe.

    Returns:
        pd.DataFrame: Merged dataframe.
    """
    merged_df = df_kl.merge(tf_family_df[['Model', 'TF family']], left_on='Target', right_on='Model', how='left')
    return merged_df

def extract_class_descriptors(df):
    """
    Extract class descriptors from the 'TF family' column.

    Args:
        df (pd.DataFrame): DataFrame containing 'TF family' column.

    Returns:
        pd.DataFrame: DataFrame with 'Class Descriptor' and new columns for class descriptors.
    """
    # Extract the class descriptors inside the braces
    df['Class Descriptor'] = df['TF family'].str.extract(r'\{(.*?)\}')

    # Split the class descriptor into parts
    descriptor_parts = df['Class Descriptor'].str.split('.')

    # Create Superclass, Class, Subclass columns
    df['Superclass'] = descriptor_parts.str[0]
    df['Class'] = descriptor_parts.str[1]
    df['Subclass'] = descriptor_parts.str[2]

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

def tag_family_by_class_descriptor(df, class_descriptors, tag_name):
    """
    Tag the dataframe based on class descriptors of varying lengths.

    Args:
        df (pd.DataFrame): DataFrame with 'Superclass', 'Class', 'Subclass' columns.
        class_descriptors (list of str): List of class descriptors to tag.
            Each class descriptor can be in 'x', 'x.y', or 'x.y.z' format.
        tag_name (str): Name of the tag column to add.

    Returns:
        pd.DataFrame: DataFrame with an additional tag column indicating family membership.
    """
    # Initialize tag column to 0
    df[tag_name] = 0

    for desc in class_descriptors:
        parts = desc.split('.')
        # Validates that parts contain only digits
        if not all(part.isdigit() for part in parts):
            print(f"Invalid class descriptor '{desc}'. It should contain only digits and periods.")
            continue  # Skip invalid descriptors
        # Build mask based on provided parts
        mask = pd.Series([True] * len(df))
        if len(parts) >= 1:
            mask &= (df['Superclass'] == parts[0])
        if len(parts) >= 2:
            mask &= (df['Class'] == parts[1])
        if len(parts) == 3:
            mask &= (df['Subclass'] == parts[2])
        df.loc[mask, tag_name] = 1
    return df

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

    def overlay_cdf(data, **kwargs):
        ax2 = plt.gca().twinx()
        sns.ecdfplot(data, ax=ax2, color='red', alpha=0.9)
        ax2.set_ylabel('CDF', color='red')
        ax2.tick_params(axis='y', labelcolor='red')
        ax2.grid(True)

    g.map(sns.histplot, 'KL distance', color='blue')
    g.map(overlay_cdf, 'KL distance')

    # Adjust titles to prevent overlapping
    g.set_titles(row_template="{row_name}", col_template="{col_name}", size=10)
    g.fig.subplots_adjust(top=0.85, hspace=0.4)

    # Use tight_layout with rect to reserve space for suptitle
    plt.tight_layout(rect=[0, 0, 1, 0.95])

    # Update the suptitle
    g.fig.suptitle(f"Histograms and CDF plots for all networks", fontsize=14)

    # Save plot
    output_path = os.path.join(output_dir, "KL_distance_distribution.png")
    plt.savefig(output_path)
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
    g.add_legend()

    # Adjust titles to prevent overlapping
    g.set_titles(row_template="{row_name}", col_template="{col_name}", size=10)
    g.fig.subplots_adjust(top=0.85, hspace=0.4)

    # Use tight_layout with rect to reserve space for suptitle
    plt.tight_layout(rect=[0, 0, 1, 0.95])

    # Update the suptitle
    g.fig.suptitle(f"Histograms by family / non-family: {family}", fontsize=14)

    # Save plot
    output_path = os.path.join(output_dir, f"KL_distribution_family_{family}.png")
    plt.savefig(output_path)
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

    # Adjust layout to prevent title overlap
    plt.tight_layout()

    # Save plot
    output_path = os.path.join(output_dir, f"KL_distance_enrichment_{family}.png")
    plt.savefig(output_path)
    plt.close()

def plot_csf_with_kstext(df_all, family, output_dir, correction_method='fdr_bh'):
    """
    Plot CDF and perform the Kolmogorov-Smirnov test for family vs non-family motifs,
    applying multiple hypothesis testing correction.

    Args:
        df_all (pd.DataFrame): DataFrame with KL distance and family tags.
        family (str): Family column name for testing.
        output_dir (str): Directory to save the output plot.
        correction_method (str): Method for multiple testing correction.
            Options include 'bonferroni', 'fdr_bh', etc.
    """
    # Get unique combinations of Dataset and Network
    groups = df_all.groupby(['Dataset', 'Network'])
    pvalues = []
    group_labels = []

    # Perform KS tests for each group and collect p-values
    for (dataset, network), group_data in groups:
        data_family = group_data[group_data[family] == 1]['KL distance'].dropna()
        data_non_family = group_data[group_data[family] == 0]['KL distance'].dropna()

        if len(data_family) == 0 or len(data_non_family) == 0:
            pvalue = np.nan
        else:
            # Perform KS test
            _, pvalue = kstest(data_family, data_non_family, alternative='less')

        pvalues.append(pvalue)
        group_labels.append((dataset, network))

    # Apply multiple testing correction
    pvalues_clean = [p for p in pvalues if not np.isnan(p)]
    if pvalues_clean:
        reject, pvals_corrected, _, _ = multipletests(pvalues_clean, method=correction_method)
    else:
        reject = []
        pvals_corrected = []
    pvals_corrected_full = []
    reject_full = []
    idx = 0
    for p in pvalues:
        if np.isnan(p):
            pvals_corrected_full.append(np.nan)
            reject_full.append(False)
        else:
            pvals_corrected_full.append(pvals_corrected[idx])
            reject_full.append(reject[idx])
            idx += 1

    # Create a DataFrame for plotting and annotations
    results_df = pd.DataFrame({
        'Dataset': [label[0] for label in group_labels],
        'Network': [label[1] for label in group_labels],
        'pvalue': pvalues,
        'pvalue_corrected': pvals_corrected_full,
        'reject': reject_full
    })

    # Plotting
    g = sns.FacetGrid(df_all, row="Dataset", col="Network", hue=family, height=4, aspect=1.5)

    def overlay_cdf(data, **kwargs):
        ax = plt.gca()
        sns.ecdfplot(data, ax=ax, **kwargs)
        ax.set_ylabel('CDF')
        ax.grid(True)

    g.map(overlay_cdf, 'KL distance')
    g.add_legend()

    # Adjust titles and layout
    g.set_titles(row_template="{row_name}", col_template="{col_name}", size=10)
    g.fig.subplots_adjust(top=0.85, hspace=0.4)

    # Use tight_layout with rect to reserve space for suptitle
    plt.tight_layout(rect=[0, 0, 1, 0.95])

    # Update the suptitle
    g.fig.suptitle(f"KS Test with Multiple Testing Correction ({correction_method})", fontsize=14)

    # Annotate plots with corrected p-values
    for ax, (_, network), pval_corr, reject_null in zip(g.axes.flatten(), group_labels, pvals_corrected_full, reject_full):
        if np.isnan(pval_corr):
            significance = 'Insufficient data'
            pval_text = 'N/A'
        else:
            significance = 'Significant' if reject_null else 'Not Significant'
            pval_text = f"{pval_corr:.3g}"
        # Place the annotation inside the plot area
        ax.text(0.05, 0.95, f"Corrected p-value: {pval_text}\n{significance}", transform=ax.transAxes,
                verticalalignment='top', fontsize=8, bbox=dict(facecolor='white', alpha=0.7, edgecolor='none'))

    # Save plot
    output_path = os.path.join(output_dir, f"KL_distribution_KS_test_{family}_corrected.png")
    plt.savefig(output_path)
    plt.close()

def perform_ks_test(df_all, family):
    """
    Perform the Kolmogorov-Smirnov test for family vs non-family motifs over the entire dataset.

    Args:
        df_all (pd.DataFrame): DataFrame with KL distance and family tags.
        family (str): Family column name for testing.

    Returns:
        dict: A dictionary containing the KS test results.
    """
    data_family = df_all[df_all[family] == 1]['KL distance'].dropna()
    data_non_family = df_all[df_all[family] == 0]['KL distance'].dropna()

    if len(data_family) == 0 or len(data_non_family) == 0:
        # Cannot perform KS test if one of the samples is empty
        result = {'ks_stat': None, 'p_value': None}
    else:
        # Perform KS test
        ks_stat, p_value = kstest(data_family, data_non_family, alternative='less')
        result = {'ks_stat': ks_stat, 'p_value': p_value}

    return result

def perform_kl_analysis(input_dir, output_dir, tf_family_file, family_files=None, family_names=None,
                        class_descriptors=None, test_all=False, test_level='superclass', correction_method='fdr_bh'):
    """
    Perform KL analysis, including plotting and family enrichment analysis.

    Args:
        input_dir (str): Directory containing KL distance data.
        output_dir (str): Directory to save all output plots.
        tf_family_file (str): Path to the TF family CSV file.
        family_files (list, optional): List of paths to family motif files.
        family_names (list, optional): List of family names.
        class_descriptors (list of str, optional): List of class descriptors to test.
            Each descriptor can be in 'x', 'x.y', or 'x.y.z' format.
        test_all (bool, optional): If True, test all superclasses or classes.
        test_level (str, optional): 'superclass' or 'class', specifies the level to test when test_all is True.
        correction_method (str): Method for multiple testing correction.
            Options include 'bonferroni', 'fdr_bh', etc.
    """
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    df_all = load_all_network_df(directory=input_dir)
    tf_family_df = load_tf_family_data(tf_family_file)

    # Merge the TF family data into the KL data
    df_all = merge_tf_family_data(df_all, tf_family_df)

    # Extract class descriptors
    df_all = extract_class_descriptors(df_all)

    family_tags = []

    if family_files and family_names:
        # Process families from files
        for family_file, family_name in zip(family_files, family_names):
            print(f"Processing family: {family_name}")
            df_all = tag_family(df_all, family_file, family_name)
            family_tags.append(f"{family_name}_tag")
    elif class_descriptors:
        # Process families specified by class descriptors
        for desc in class_descriptors:
            if not desc.replace('.', '').isdigit():
                print(f"Invalid class descriptor '{desc}'. It should contain only digits and periods.")
                continue
            tag_name = f"class_{desc.replace('.', '_')}_tag"
            print(f"Processing class descriptor: {desc}")
            df_all = tag_family_by_class_descriptor(df_all, [desc], tag_name)
            family_tags.append(tag_name)
    elif test_all:
        # Get all unique levels in the data based on test_level
        if test_level == 'superclass':
            levels = df_all['Superclass'].dropna().unique()
            for level in levels:
                tag_name = f"superclass_{level}_tag"
                print(f"Processing superclass: {level}")
                # Tag based on Superclass
                df_all[tag_name] = (df_all['Superclass'] == level).astype(int)
                family_tags.append(tag_name)
        elif test_level == 'class':
            levels = df_all[['Superclass', 'Class']].dropna().drop_duplicates()
            for _, row in levels.iterrows():
                level = f"{row['Superclass']}.{row['Class']}"
                tag_name = f"class_{level.replace('.', '_')}_tag"
                print(f"Processing class: {level}")
                df_all = tag_family_by_class_descriptor(df_all, [level], tag_name)
                family_tags.append(tag_name)
        else:
            raise ValueError("Invalid test_level. Must be 'superclass' or 'class'.")
    else:
        print("No families or class descriptors specified for analysis.")
        return

    if test_all:
        # For test_all case, collect KS test results
        ks_results = []

        for tag in family_tags:
            print(f"Performing KS test for {tag}")
            ks_result = perform_ks_test(df_all, family=tag)
            ks_result['Family'] = tag
            ks_results.append(ks_result)

        # Apply multiple testing correction
        p_values = [result['p_value'] for result in ks_results if result['p_value'] is not None]
        families_with_pvalues = [result['Family'] for result in ks_results if result['p_value'] is not None]

        if p_values:
            reject, pvals_corrected, _, _ = multipletests(p_values, method=correction_method)
            # Update ks_results with corrected p-values and rejection decisions
            idx = 0
            for result in ks_results:
                if result['p_value'] is not None:
                    result['p_value_corrected'] = pvals_corrected[idx]
                    result['reject'] = reject[idx]
                    idx += 1
                else:
                    result['p_value_corrected'] = None
                    result['reject'] = False
        else:
            for result in ks_results:
                result['p_value_corrected'] = None
                result['reject'] = False

        # Save ks_results to CSV
        ks_results_df = pd.DataFrame(ks_results)
        output_csv = os.path.join(output_dir, f'ks_test_results_{test_level}.csv')
        ks_results_df.to_csv(output_csv, index=False)
        print(f"KS test results saved to {output_csv}")
    else:
        for tag in family_tags:
            print(f"Processing {tag}")
            df_all_tagged = df_all.copy()
            # Proceed to plotting and analysis with df_all_tagged
            plot_kl_distributions(df_all_tagged, output_dir)
            plot_kl_distributions_with_motif_family(df_all_tagged, family=tag, output_dir=output_dir)
            sorted_df = df_all_tagged.sort_values(by=["KL distance"])
            enrichment_df = get_enrichment_df(sorted_df, tag)
            plot_enrichment(enrichment_df, family=tag, output_dir=output_dir)
            plot_csf_with_kstext(df_all_tagged, family=tag, output_dir=output_dir, correction_method=correction_method)
