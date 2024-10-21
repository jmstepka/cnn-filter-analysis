from Bio import SeqIO
from Bio.Seq import Seq
from Bio import motifs
import numpy as np
import pandas as pd
import os
import random


def read_sequences(file, filetype="fasta"):
    """
    Function to read sequences from a FASTA file.
    
    :param file: File with sequences
    :param filetype: File type, default is fasta
    :return: List of sequences
    """
    sequences = []
    for seq_record in SeqIO.parse(file, filetype):
        sequences.append(seq_record.seq.upper())
    return sequences


def one_hot_encode(seq):
    """
    Function to one-hot encode a nucleotide sequence.
    
    :param seq: Sequence to encode
    :return: One-hot encoded matrix
    """
    if "N" in seq:
        new_seq = ""
        for letter in seq:
            if letter == "N":
                new_seq += random.choice("ACGT")
            else:
                new_seq += letter
        seq = new_seq
    
    mapping = dict(zip("ACGT", range(4)))
    seq2 = [mapping[i] for i in seq]
    return np.eye(4)[seq2]


def calculate_filter_mean(filter_matrix):
    """
    Function to calculate the mean value of filter components.
    
    :param filter_matrix: Filter matrix
    :return: Mean value of the filter matrix elements
    """
    return np.mean(filter_matrix)


def scan_sequences(sequences_file, filters_file, threshold=1e-5):
    """
    Function to scan sequences with a filter. Returns a dictionary where the key
    is the filter name and the value is a list of the best subsequences for each sequence.
    
    :param sequences_file: File with sequences from the training set
    :param filters_file: File with filter values
    :param threshold: Threshold for the mean filter value, above which the filter is considered
    :return: Dictionary with the best subsequences for each filter
    """
    
    best_sequences = {}
    sequences = read_sequences(sequences_file)
    
    nucleotides = {0: "A", 1: "C", 2: "G", 3: "T"}

    for sequence in sequences:
        matrix = one_hot_encode(sequence).transpose()
        
        with open(filters_file, "r") as file:
            lines = file.readlines()
            f = 0
            for i in range(1, len(lines)-1, 5):
                filter_matrix = np.loadtxt(lines[i:i+4])
                filter_mean = calculate_filter_mean(filter_matrix)
                
                # Filter filters based on the mean value
                if filter_mean > threshold:
                    best = float('-inf')
                    for j in range(np.shape(matrix)[1]):
                        if np.shape(matrix[:, j:j+np.shape(filter_matrix)[1]])[1] == np.shape(filter_matrix)[1]:
                            submatrix = matrix[:, j:j+np.shape(filter_matrix)[1]]
                            convolution = filter_matrix * submatrix
                            total = np.sum(convolution)
                            if total > best:
                                best = total
                                best_seq = submatrix
                                where = np.argwhere(best_seq)
                                seq = "*" * np.shape(best_seq)[1]
                                for index in where:
                                    seq = seq[:index[1]] + nucleotides[index[0]] + seq[index[1]+1:]
                    if f"filter_{f}" in best_sequences.keys():
                        tmp = best_sequences[f"filter_{f}"]
                        tmp.append((seq, float(best)))
                    else:
                        best_sequences[f"filter_{f}"] = [(seq, float(best))]
                f += 1
    return best_sequences


def choose_best_subseq(dictionary, n=100):
    """
    Function to choose the n best subsequences for each filter.
    
    :param dictionary: Dictionary with subsequences and their convolution scores
    :param n: Number of best subsequences to choose
    :return: Dictionary with the best n subsequences
    """
    dictionary_new = {}
    for key in dictionary:
        sorted_list = sorted(dictionary[key], key=lambda tup: tup[1], reverse=True)
        tmp = []
        for i in range(n):
            tmp.append(sorted_list[i][0])
        dictionary_new[key] = tmp

    return dictionary_new


def calculate_pfm_matrices(dictionary, dataset, network, class_type, output_dir):
    """
    Function to create a PFM matrix based on a dictionary where
    the key is the filter name and the value is a list of the best
    n subsequences. Saves the created PFM matrices for each filter
    to a file in the specified output location.

    :param dictionary: Dictionary where the key is the filter name 
                       and the value is a list of the best n subsequences
    :param dataset: Dataset name
    :param network: Network type
    :param class_type: Class type
    :param output_dir: Output directory where the PFM matrices should be saved
    """
    
    # Path to the output file
    output_file_path = os.path.join(output_dir, f"{dataset}_{network}_{class_type}_pfm.txt")
    
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    with open(output_file_path, "w") as fileout:
        fileout.write("MEME version 4\n\nALPHABET= ACGT\n\nstrands: + -\n\nBackground letter frequencies\nA 0.25 C 0.25 G 0.25 T 0.25\n\n")
        # For each filter
        for filter_name in dictionary:
            fileout.write(f"MOTIF {str(filter_name)}\nletter-probability matrix: alength= 4 w= 19\n")
            instances = []
            # Creating Seq instances for each subsequence 
            for seq in dictionary[filter_name]:
                instances.append(Seq(seq))
            # Creating PFM matrix for the set of subsequences for each filter
            m = motifs.create(instances)
            pfm_matrix = m.counts
            # Converting Motif instances to a numpy matrix
            pfm_list = []
            for nucleotide in ["A", "C", "G", "T"]:
                pfm_list.append(list(pfm_matrix[nucleotide]))
            pfm_list = np.matrix(pfm_list).transpose()
            np.savetxt(fileout, pfm_list)
            fileout.write("\n")
    
    print(f"Saved PFM matrix to file: {output_file_path}")


def create_pfm_matrices(input_dir, output_dir, threshold=1e-5):
    """
    Function to generate PFM matrices based on files from the `filters` and `training_data` folders.
    
    :param input_dir: Input location where the `filters` and `training_data` folders are located.
    :param output_dir: Output location where the PFM matrices should be saved.
    :param threshold: Threshold for the mean filter value, above which the filter is considered.
    """
    
    filters_dir = os.path.join(input_dir, "filters")
    training_data_dir = os.path.join(input_dir, "training_data")

    # Iterate over files in the filters folder
    for filters_file in os.listdir(filters_dir):
        if filters_file.endswith("_filter.txt"):
            # Extract information from the file name
            dataset, network_type = filters_file.replace("_filter.txt", "").split("_")
            filters_file_path = os.path.join(filters_dir, filters_file)

            # Iterate over files in the training_data folder
            for training_data_file in os.listdir(training_data_dir):
                if training_data_file.endswith(".fa"):
                    # Extract information from the file name
                    dataset_td, class_type = training_data_file.replace(".fa", "").split("_")

                    # Ensure that the dataset from filters matches the dataset from training_data
                    if dataset == dataset_td:
                        training_data_file_path = os.path.join(training_data_dir, training_data_file)
                        
                        # Scan sequences and generate PFM
                        best_sequences = scan_sequences(training_data_file_path, filters_file_path, threshold=threshold)
                        best_sequences_n = choose_best_subseq(best_sequences)
                        calculate_pfm_matrices(best_sequences_n, dataset, network_type, class_type, output_dir)

    print("PFM matrix creation process completed.")
