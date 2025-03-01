import numpy as np
import pandas as pd
import os

def add_pseudocount(freq, sites, pseudocount=0.1):
    """
    Adds pseudocount to frequency assuming uniform background.

    Args:
        freq (float): Frequency value to adjust.
        sites (int): Number of binding sites.
        pseudocount (float): Pseudocount to add for adjustment.

    Returns:
        float: Adjusted frequency with pseudocount.
    """
    return ((freq * sites) + (pseudocount * 0.25)) / (sites + pseudocount)

def calculate_kl(query_freq, target_freq, epsilon=1e-10):
    """
    Calculate the Kullback-Leibler (KL) divergence between two frequencies.

    Args:
        query_freq (float): Query frequency value.
        target_freq (float): Target frequency value.
        epsilon (float): Small value to avoid log(0).

    Returns:
        float: Average KL divergence between the two frequencies.
    """
    query_freq = max(query_freq, epsilon)  # Prevent log(0)
    target_freq = max(target_freq, epsilon)
    
    avg_kull = ((query_freq * np.log10(query_freq / target_freq)) + 
                (target_freq * np.log10(target_freq / query_freq))) / 2
    return avg_kull

def filter_motifs_by_ic(ic_csv_file, ic_threshold):
    """
    Filters motifs based on the Information Content (IC) threshold.

    Args:
        ic_csv_file (str): Path to the IC CSV file.
        ic_threshold (float): IC threshold value.

    Returns:
        list: List of motifs that meet the IC threshold.
    """
    ic_df = pd.read_csv(ic_csv_file)
    filtered_motifs = ic_df[ic_df["IC"] >= ic_threshold]["Filter"].tolist()
    return filtered_motifs

def load_pfm_matrices(pfm_folder, dataset, network, class_type):
    """
    Loads PFM matrices for a given dataset, network, and class type.

    Args:
        pfm_folder (str): Folder path containing the PFM matrices.
        dataset (str): Dataset name.
        network (str): Network name.
        class_type (str): Class type.

    Returns:
        dict: Dictionary of PFM matrices with motif names as keys.
    """
    pfm_file_path = os.path.join(pfm_folder, f"{dataset}_{network}_{class_type}_pfm.txt")
    pfm_matrices = {}
    
    if os.path.exists(pfm_file_path):
        with open(pfm_file_path) as pfm_file:
            motifs = pfm_file.read().split("\nMOTIF ")[1:]
            for motif in motifs:
                qname = motif.split("\n")[0]
                matrix = np.loadtxt(motif.split("\n")[2:])
                pfm_matrices[qname] = matrix
    return pfm_matrices

def create_ppm_matrix(pfm_matrix):
    """
    Converts a PFM matrix into a PPM matrix.

    Args:
        pfm_matrix (np.array): Position Frequency Matrix (PFM).

    Returns:
        np.array: Position Probability Matrix (PPM).
    """
    ppm_matrix = np.empty([pfm_matrix.shape[0], pfm_matrix.shape[1]])
    for (x, y), value in np.ndenumerate(pfm_matrix):
        total = np.sum(pfm_matrix[x])
        ppm_matrix[x][y] = value / total
    return ppm_matrix

def extract_kmer_from_ppm(ppm_matrix, start, end):
    """
    Extracts a k-mer (submatrix) from the PPM matrix.

    Args:
        ppm_matrix (np.array): Position Probability Matrix (PPM).
        start (int): Start position of the k-mer.
        end (int): End position of the k-mer.

    Returns:
        np.array: Submatrix representing the k-mer.
    """
    return ppm_matrix[start:end]

def calculate_kl_distance(ic_csv_folder, pfm_folder, output_folder, hocomoco_models_path, ic_threshold=6.5):
    """
    Calculates KL divergence for motifs filtered by IC threshold and saves the results.

    Args:
        ic_csv_folder (str): Folder containing IC CSV files.
        pfm_folder (str): Folder containing PFM matrices.
        output_folder (str): Folder to save the KL divergence results.
        hocomoco_models_path (str): Path to HOCOMOCO models in MEME format.
        ic_threshold (float): IC threshold for motif filtering.

    Returns:
        None: Saves KL divergence results to CSV files in the output folder.
    """
    alphabet_size = 4

    for ic_csv_file in os.listdir(ic_csv_folder):
        if ic_csv_file.endswith(".csv"):
            filename_parts = ic_csv_file.replace(".csv", "").split("_")
            dataset, network, class_type = filename_parts

            print(f"Processing {dataset}, {network}, {class_type}")
            ic_csv_file_path = os.path.join(ic_csv_folder, ic_csv_file)

            # Filter motifs by IC threshold using the CSV file
            filter_names = filter_motifs_by_ic(ic_csv_file_path, ic_threshold)

            pfm_matrices = load_pfm_matrices(pfm_folder, dataset, network, class_type)

            results_dict = {
                "Dataset": [], "Network": [], "Class": [], "Query": [],
                "Target": [], "Offset": [], "KL distance": []
            }

            # Iterate over each filtered motif name from the CSV file
            for filter_name in filter_names:
                filter_parts = filter_name.split("_")
                filter_num = filter_parts[4]
                start_pos = int(filter_parts[5])
                end_pos = int(filter_parts[6])

                # Retrieve the corresponding PFM matrix for the filter
                pfm_matrix = pfm_matrices.get(f"filter_{filter_num}")
                if pfm_matrix is None:
                    continue 

                # Convert PFM to PPM matrix
                ppm_matrix = create_ppm_matrix(pfm_matrix)

                # Extract the kmer from the PPM matrix
                kmer_ppm = extract_kmer_from_ppm(ppm_matrix, start_pos, end_pos)

                length = kmer_ppm.shape[0]  # Length of the kmer

                # Iterate over target motifs in HOCOMOCO
                with open(hocomoco_models_path) as tfile:
                    tmotifs = tfile.read().split("\nMOTIF ")[1:]
                    for tmotif in tmotifs:
                        tname = tmotif.split("\n")[0]
                        tmatrix = np.loadtxt(tmotif.split("\nURL")[0].split("\n")[2:])
                        nsites = int(tmotif.split("\n")[1].split("nsites= ")[-1])
                        tlength = tmatrix.shape[0]
                        if tlength < length:
                            continue

                        rc_tmatrix = np.flip(tmatrix, axis=0)  # Reverse complement of the target matrix

                        distances_dict = {}

                        # Sliding window on the target motif
                        for iq in range(0, length):
                            for it in range(iq, iq + tlength - length + 1):
                                offset = it - iq                        
                                distance = 0
                                rc_distance = 0

                                for ia in range(0, alphabet_size):                                                            
                                    q_freq = kmer_ppm[iq, ia]
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
                            results_dict["Query"].append(filter_name)
                            results_dict["Target"].append(tname)
                            results_dict["Offset"].append(offset)
                            results_dict["KL distance"].append(lowest_scores_dict[offset])


            if not os.path.exists(output_folder):
                os.makedirs(output_folder)
            output_file = os.path.join(output_folder, f"{dataset}_{network}_{class_type}_KLdist.csv")
            results = pd.DataFrame(results_dict)
            results.to_csv(output_file, sep="\t", index=False)