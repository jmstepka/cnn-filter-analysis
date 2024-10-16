import numpy as np
import pandas as pd
import os

def add_pseudocount(freq, sites, pseudocount=0.1):
    """
    Add pseudocount assuming uniform background.
    """
    return ((freq * sites) + (pseudocount * 0.25)) / (sites + pseudocount)

def calculate_kl(query_freq, target_freq):
    """
    Calculate Kullback-Leibler divergence between two positions.
    """
    avg_kull = ((query_freq * np.log10(query_freq / target_freq)) + 
                (target_freq * np.log10(target_freq / query_freq))) / 2
    return avg_kull

def filter_motifs_by_ic(ic_csv_file, ic_threshold):
    """
    Filter motifs based on Information Content (IC) threshold from a given CSV file.
    """
    ic_df = pd.read_csv(ic_csv_file)
    filtered_motifs = ic_df[ic_df["IC"] >= ic_threshold]["Filter"].tolist()
    return filtered_motifs

def load_pfm_matrices(pfm_folder, filter_names):
    """
    Load PFM matrices for the given filter names from the PFM files.
    """
    pfm_matrices = {}
    
    for filter_name in filter_names:
        pfm_file_path = os.path.join(pfm_folder, f"{filter_name}.txt")
        if os.path.exists(pfm_file_path):
            with open(pfm_file_path) as pfm_file:
                motifs = pfm_file.read().split("\nMOTIF ")[1:]
                for motif in motifs:
                    qname = motif.split("\n")[0]
                    matrix = np.loadtxt(motif.split("\n")[2:])
                    pfm_matrices[qname] = matrix
    return pfm_matrices

def calculate_kl_distance(ic_csv_folder, pfm_folder, output_folder, hocomoco_models_path, ic_threshold=6.5):
    """
    Calculate KL divergence for filtered motifs based on IC threshold from CSV files.
    """
    alphabet_size = 4
    pseudocount = 0.1

    # Iterate through all CSV files in the IC folder
    for ic_csv_file in os.listdir(ic_csv_folder):
        if ic_csv_file.endswith(".csv"):
            # Extract dataset, network, and class from the filename
            filename_parts = ic_csv_file.replace(".csv", "").split("_")
            dataset, network, class_type = filename_parts

            # Define paths and print progress
            print(f"Processing {dataset}, {network}, {class_type}")
            ic_csv_file_path = os.path.join(ic_csv_folder, ic_csv_file)
            pfm_input_path = os.path.join(pfm_folder, f"{dataset}_{network}_{class_type}_pfm.txt")

            # Filter motifs by IC threshold using the CSV file
            filter_names = filter_motifs_by_ic(ic_csv_file_path, ic_threshold)

            # Load corresponding PFM matrices
            pfm_matrices = load_pfm_matrices(pfm_input_path, filter_names)

            results_dict = {
                "Dataset": [], "Network": [], "Class": [], "Query": [],
                "Target": [], "Offset": [], "KL distance": []
            }

            # Iterate over each query motif (PFM matrix)
            for qname, matrix in pfm_matrices.items():
                length = len(matrix)

                # Iterate over target motifs in HOCOMOCO
                with open(hocomoco_models_path) as tfile:
                    tmotifs = tfile.read().split("\nMOTIF ")[1:]
                    for tmotif in tmotifs:
                        tname = tmotif.split("\n")[0]
                        tmatrix = np.loadtxt(tmotif.split("\nURL")[0].split("\n")[2:])
                        nsites = int(tmotif.split("\n")[1].split("nsites= ")[-1])
                        tlength = len(tmatrix)

                        rc_tmatrix = np.flip(tmatrix)

                        distances_dict = {}

                        # Sliding window on the target motif
                        for iq in range(0, length):
                            for it in range(iq, iq + tlength - length + 1):
                                offset = it - iq                        
                                distance = 0
                                rc_distance = 0

                                for ia in range(0, alphabet_size):                                                
                                    q_freq = matrix[iq, ia]
                                    t_freq = tmatrix[it, ia]
                                    rc_t_freq = rc_tmatrix[it, ia]

                                    t_freq = add_pseudocount(t_freq, nsites)
                                    rc_t_freq = add_pseudocount(rc_t_freq, nsites)

                                    avg_kull = calculate_kl(q_freq, t_freq)
                                    avg_kull_rc = calculate_kl(q_freq, rc_t_freq)

                                    distance += avg_kull
                                    rc_distance += avg_kull_rc

                                distance = -1 * distance
                                rc_distance = -1 * rc_distance

                                if offset in distances_dict:
                                    distances_dict[offset] += distance
                                    distances_dict[-(offset+1)] += rc_distance
                                else:
                                    distances_dict[offset] = distance
                                    distances_dict[-(offset+1)] = rc_distance

                        offset_dist_list = [[o, d] for o, d in distances_dict.items()]
                        offset_dist_list_sorted = sorted(offset_dist_list, key=lambda x: x[1], reverse=True)
                        lowest_score = offset_dist_list_sorted[0][1]
                        lowest_scores_dict = {offset_dist_list_sorted[0][0]: offset_dist_list_sorted[0][1]}

                        for o, d in offset_dist_list_sorted:
                            if d == lowest_score:
                                lowest_scores_dict[o] = d
                            else:
                                break

                        for offset in lowest_scores_dict.keys():
                            results_dict["Dataset"].append(dataset)
                            results_dict["Network"].append(network)
                            results_dict["Class"].append(class_type)
                            results_dict["Query"].append(qname)
                            results_dict["Target"].append(tname)
                            results_dict["Offset"].append(offset)
                            results_dict["KL distance"].append(lowest_scores_dict[offset])

            # Save results to CSV file
            if not os.path.exists(output_folder):
                os.makedirs(output_folder)
            output_file = os.path.join(output_folder, f"{dataset}_{network}_{class_type}_KLdist.csv")
            results = pd.DataFrame(results_dict)
            results.to_csv(output_file, sep="\t", index=False)
