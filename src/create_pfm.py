from Bio import SeqIO
from Bio.Seq import Seq
from Bio import motifs
import numpy as np
import pandas as pd
import os
import random


def read_sequences(file, filetype="fasta"):
    """
    Funkcja wczytująca sekwencje z pliku FASTA.
    
    :param file: Plik z sekwencjami
    :param filetype: Typ pliku, domyślnie fasta
    :return: Lista sekwencji
    """
    sequences = []
    for seq_record in SeqIO.parse(file, filetype):
        sequences.append(seq_record.seq.upper())
    return sequences


def one_hot_encode(seq):
    """
    Funkcja kodująca sekwencję nukleotydów w sposób one-hot encoding.
    
    :param seq: Sekwencja do zakodowania
    :return: Macierz one-hot encoding
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
    Funkcja obliczająca średnią wartość składowych filtra.
    
    :param filter_matrix: Macierz filtra
    :return: Średnia wartość elementów macierzy filtra
    """
    return np.mean(filter_matrix)


def scan_sequences(sequences_file, filters_file, threshold=1e-5):
    """
    Funkcja skanująca filtrem sekwencje. Zwraca słownik, w którym kluczem
    jest nazwa filtra, a wartością lista najlepszych podsekwencji dla każdej sekwencji.
    
    :param sequences_file: Plik z sekwencjami ze zbioru treningowego
    :param filters_file: Plik z wartościami filtra
    :param threshold: Próg średniej wartości filtra, powyżej którego filtr jest brany pod uwagę
    :return: Słownik z najlepszymi podsekwencjami dla każdego filtra
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
                
                # Filtruj filtry na podstawie średniej wartości
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
    Funkcja wybierająca n najlepszych podsekwencji dla każdego filtra.
    
    :param dictionary: Słownik z podsekwencjami i ich wynikami splotu
    :param n: Liczba najlepszych podsekwencji do wybrania
    :return: Słownik z najlepszymi n podsekwencjami
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
    Funkcja tworzy macierz PFM na podstawie słownika, w którym
    kluczem jest nazwa filtra, a wartością lista najlepszych
    n podsekwencji. Zapisuje utworzone macierze PFM dla każdego
    filtra do pliku w podanej lokalizacji wyjściowej.

    :param dictionary: Słownik, w którym kluczem jest nazwa filtra, 
                       a wartością lista najlepszych n podsekwencji
    :param dataset: Nazwa datasetu
    :param network: Typ sieci (network type)
    :param class_type: Typ klasy
    :param output_dir: Katalog wyjściowy, w którym mają być zapisane macierze PFM
    """
    
    # Ścieżka do pliku wyjściowego
    output_file_path = os.path.join(output_dir, f"{dataset}_{network}_{class_type}_pfm.txt")
    
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    with open(output_file_path, "w") as fileout:
        fileout.write("MEME version 4\n\nALPHABET= ACGT\n\nstrands: + -\n\nBackground letter frequencies\nA 0.25 C 0.25 G 0.25 T 0.25\n\n")
        # Dla każdego filtra
        for filter_name in dictionary:
            fileout.write(f"MOTIF {str(filter_name)}\nletter-probability matrix: alength= 4 w= 19\n")
            instances = []
            # Tworzenie instancji Seq dla każdej podsekwencji 
            for seq in dictionary[filter_name]:
                instances.append(Seq(seq))
            # Tworzenie macierzy PFM dla zbioru podsekwencji dla każdego filtra
            m = motifs.create(instances)
            pfm_matrix = m.counts
            # Zamiana instancji Motif na macierz numpy
            pfm_list = []
            for nucleotide in ["A", "C", "G", "T"]:
                pfm_list.append(list(pfm_matrix[nucleotide]))
            pfm_list = np.matrix(pfm_list).transpose()
            np.savetxt(fileout, pfm_list)
            fileout.write("\n")
    
    print(f"Zapisano macierz PFM do pliku: {output_file_path}")


def create_pfm_matrices(input_dir, output_dir, threshold=1e-5):
    """
    Funkcja generująca macierze PFM na podstawie plików z folderów `filters` i `training_data`.
    
    :param input_dir: Lokalizacja wejściowa, w której znajdują się foldery `filters` i `training_data`.
    :param output_dir: Lokalizacja wyjściowa, gdzie mają być zapisane macierze PFM.
    :param threshold: Próg średniej wartości filtra, powyżej którego filtr jest brany pod uwagę.
    """
    
    filters_dir = os.path.join(input_dir, "filters")
    training_data_dir = os.path.join(input_dir, "training_data")

    # Iteracja po plikach w folderze filters
    for filters_file in os.listdir(filters_dir):
        if filters_file.endswith("_filter.txt"):
            # Wyciągnięcie informacji z nazwy pliku
            dataset, network_type = filters_file.replace("_filter.txt", "").split("_")
            filters_file_path = os.path.join(filters_dir, filters_file)

            # Iteracja po plikach w folderze training_data
            for training_data_file in os.listdir(training_data_dir):
                if training_data_file.endswith(".fa"):
                    # Wyciągnięcie informacji z nazwy pliku
                    dataset_td, class_type = training_data_file.replace(".fa", "").split("_")

                    # Upewnienie się, że dataset z filters pasuje do dataset z training_data
                    if dataset == dataset_td:
                        training_data_file_path = os.path.join(training_data_dir, training_data_file)
                        
                        # Skany sekwencji i generowanie PFM
                        best_sequences = scan_sequences(training_data_file_path, filters_file_path, threshold=threshold)
                        best_sequences_n = choose_best_subseq(best_sequences)
                        calculate_pfm_matrices(best_sequences_n, dataset, network_type, class_type, output_dir)

    print("Proces tworzenia macierzy PFM zakończony.")
