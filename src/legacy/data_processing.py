
import pandas as pd
import numpy as np
import logomaker as lm

from pathlib import Path

species = ["hg38_5k", "pt06_5k", "rm10_5k"]
variants = ["last", "again_last"]
classes = ["pa", "pi", "na", "ni"]
common_dataset_names = ["promoter_active",
                        "promoter_inactive",
                        "nonpromoter_active",
                        "nonpromoter_inactive"]


def get_network_dataset_dict():
    network_dataset_dict = {}
    for s in species:
        for n in variants:
            network_name = f"{s}_{n}"
            network_dataset_dict[network_name] = []
            for d in classes:
                network_dataset_dict[network_name].append(f"{s}_{d}")

    return network_dataset_dict


def load_one_network_df(network, directory="output"):
    dfs = []
    for c in classes:
        dfs.append(pd.read_csv(Path(directory) / f"{network}_{c}_KLdist.csv", sep="\t"))
    return pd.concat(dfs)


def load_all_network_df(directory="output_50"):
    dfs = []
    for network in get_network_dataset_dict():
        for c in classes:
            dfs.append(pd.read_csv(Path(directory) / f"{network}_{c}_KLdist.csv", sep="\t"))
    df = pd.concat(dfs)
    df["network_type"] = df["Network"].str[8:]
    df["KL distance"] = -1 * df["KL distance"]
    return df


def create_ppm_matrix(pfm_matrix):
    """
    by J. Smolik
    Funkcja tworzaca macierz ppm na podstawie macierzy pfm.

    :param pfm_matrix: macierz pfm
    """
    
    ppm_matrix = np.empty([pfm_matrix.shape[0], pfm_matrix.shape[1]])
    for (x,y), value in np.ndenumerate(pfm_matrix):
        total = np.sum(pfm_matrix[x])
        ppm_matrix[x][y] = value/total
    return ppm_matrix


def get_matrix(row):
    pfm_file = "pfm_filters/{}_5k_{}/{}_pfm.txt".format(row["Species"], row["Dataset"], row["Network"])
    filter = "_".join(row["Query"].split("_")[:2])
    kmer_start, kmer_end = int(row["Query"].split("_")[2]), int(row["Query"].split("_")[3])
    with open(pfm_file) as file:
        lines = file.readlines()
        for i, l in enumerate(lines):
            if l == f"MOTIF {filter}\n":
                print(pfm_file, filter, i)
                pfm = np.loadtxt(lines[i+2+kmer_start:i+2+kmer_end])
                break
    ppm_matrix = create_ppm_matrix(pfm)
    ppm = pd.DataFrame(ppm_matrix, columns = ["A", "C", "G", "T"])
    info_mat = lm.transform_matrix(ppm, 
                    to_type='information', 
                    from_type='probability')  
    return info_mat


def get_tf_matrix(row):
    filepath = "data/HOCOMOCOv11_core_HUMAN_mono_meme_format.meme"
    motif = row["Target"]
    offset = row["Offset"]
    print(offset)
    reverse = False
    if offset < 0:
        reverse = True
        offset = -offset - 1
    start=0
    with open(filepath) as file:
        lines = file.readlines()
        for i, l in enumerate(lines):
            if l == f"MOTIF {motif}\n":
                start = i
                
            if l.split(' ')[0] == "URL" and start != 0:
                end = i
                break
    pfm = np.loadtxt(lines[start+2:end])
    if reverse:
        pfm = np.flip(pfm)
    pfm = pfm[offset:offset+5]
    ppm_matrix = create_ppm_matrix(pfm)
    ppm = pd.DataFrame(ppm_matrix, columns = ["A", "C", "G", "T"])
    info_mat = lm.transform_matrix(ppm, 
                    to_type='information', 
                    from_type='probability')  
    return info_mat


def tag_family(df, path, name):
    family_tfs = pd.read_csv(f"data/{path}", sep = "\t")
    family_tfs[f"{name}_tag"] = 1
    family_df = family_tfs[["Model", f"{name}_tag"]]
    tagged = df.merge(family_df, left_on="Target", right_on="Model", how="left")
    tagged[f"{name}_tag"].fillna(0, inplace=True)
    return tagged


def get_enrichment_df(sorted_df, tag_name, shortening_thresh=1.5):
    dfs = []
    for network in get_network_dataset_dict():
        enrichment_data = []
        dist_data = []
        total = 0
        stripe_count = 0
        for r in sorted_df[sorted_df["Network"] == network].iterrows():
            r = r[1]
            stripe = r[tag_name]
            total += 1
            stripe_count += stripe
            enrichment_data.append(stripe_count/total)
            dist_data.append(r["KL distance"])
        dict = {"enrichment": enrichment_data, "KL distance": dist_data}
        df = pd.DataFrame(dict)
        df["network"] = network
        dfs.append(df)
    enrichment_df = pd.concat(dfs, axis=0)
    enrichment_df_shortened = pd.concat([enrichment_df[enrichment_df["KL distance"] <= shortening_thresh],
                                    enrichment_df[enrichment_df["KL distance"] > shortening_thresh][::20]])
    return enrichment_df_shortened