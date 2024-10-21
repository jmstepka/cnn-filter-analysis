import os
import numpy as np
import csv

def create_pseudo_counts_matrix(pfm_matrix, epsilon=1):
    """
    Create a pseudo counts matrix from a PFM matrix.
    
    Args:
        pfm_matrix (numpy.ndarray): Position Frequency Matrix.
        epsilon (float): Small constant to avoid zero counts.
    
    Returns:
        numpy.ndarray: Pseudo counts matrix.
    """
    pseudo_counts_matrix = np.empty([pfm_matrix.shape[0], pfm_matrix.shape[1]])
    for (x, y), value in np.ndenumerate(pfm_matrix):
        pseudo_counts_matrix[x][y] = (value + (epsilon / 4)) / (100 + epsilon)
    return pseudo_counts_matrix

def create_ppm_matrix(pseudo_counts_matrix):
    """
    Create a PPM matrix from a pseudo counts matrix.
    
    Args:
        pseudo_counts_matrix (numpy.ndarray): Pseudo counts matrix.
    
    Returns:
        numpy.ndarray: Position Probability Matrix.
    """
    ppm_matrix = np.empty([pseudo_counts_matrix.shape[0], pseudo_counts_matrix.shape[1]])
    for (x, y), value in np.ndenumerate(pseudo_counts_matrix):
        total = np.sum(pseudo_counts_matrix[x])
        ppm_matrix[x][y] = value / total
    return ppm_matrix

def ic_content(ppm_matrix):
    """
    Calculate the Information Content (IC) from a PPM matrix.
    
    Args:
        ppm_matrix (numpy.ndarray): Position Probability Matrix.
    
    Returns:
        float: Total Information Content.
    """
    ic_matrix = np.empty([ppm_matrix.shape[0], ppm_matrix.shape[1]])
    for r in range(len(ppm_matrix)):
        uncertainty = 0
        for value in ppm_matrix[r]:
            if value > 0:
                uncertainty += value * np.log2(value)
        uncertainty = (-1) * uncertainty
        ic_final = 2 - uncertainty
        row = []
        for value in ppm_matrix[r]:
            row.append(value * ic_final)
        ic_matrix[r] = row
    total_ic = np.sum(ic_matrix)
    return total_ic

def calculate_ic(input_folder, output_folder, k=7):
    """
    Calculate the Information Content (IC) for each motif in the input folder and save the results to the output folder.
    
    Args:
        input_folder (str): Path to the folder containing input PFM files.
        output_folder (str): Path to the folder where output CSV files will be saved.
        k (int): Length of the k-mer to be analyzed.
    """
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    for filename in os.listdir(input_folder):
        if filename.endswith("_pfm.txt"):
            filepath = os.path.join(input_folder, filename)
            dataset, network, class_type = filename.split("_")[:3]

            output_file = os.path.join(output_folder, f"{dataset}_{network}_{class_type}.csv")
            with open(filepath) as pfm_file:
                lines_tmp = pfm_file.readlines()

                # Skip the first 9 lines
                lines = lines_tmp[9:]

                motifs = []
                current_motif = []
                current_filter_num = None
                for line in lines:
                    if line.startswith("MOTIF"):
                        # Extract filter number
                        current_filter_num = line.split()[1].split("_")[1]
                        if current_motif:
                            motifs.append((current_filter_num, current_motif))
                            current_motif = []
                    elif line.strip() and "letter-probability matrix:" not in line:
                        current_motif.append(line.strip())
                if current_motif:
                    motifs.append((current_filter_num, current_motif))  # add the last motif

                with open(output_file, "w", newline="") as csv_output:
                    writer_ic = csv.writer(csv_output)
                    writer_ic.writerow(["Filter", "IC"])

                    # Iterate through each motif
                    for filter_num, motif in motifs:
                        weights_matrix = np.loadtxt(motif)
                        num_rows, _ = weights_matrix.shape

                        for i in range(0, num_rows - k + 1):
                            kmer_freq = weights_matrix[i:(i + k)]
                            kmer_pseudo_counts = create_pseudo_counts_matrix(kmer_freq)
                            kmer_ppm = create_ppm_matrix(kmer_pseudo_counts)
                            info_content = ic_content(kmer_ppm)

                            # Save filter_{filter_num} to CSV file
                            writer_ic.writerow([f"{dataset}_{network}_{class_type}_filter_{filter_num}_{i}_{i + k}", info_content])

