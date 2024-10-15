from Bio import SeqIO
from Bio.Seq import Seq
from Bio import motifs
import numpy as np
import pandas as pd
import os
import random


def read_sequences(file, filetype = "fasta"):
    """
    Funkcja wczytujaca sekwencje z pliku
    
    :param file: plik z sekwencjami
    :param filetype: typ pliku, domyslnie fasta
    """
    
    sequences = []
    
    for seq_record in SeqIO.parse(file, filetype):
        sequences.append(seq_record.seq.upper())
    
    return sequences


def one_hot_encode(seq):
    """
    Funkcja kodujaca sekwencje nukleotydow w sposob one hot encoding
    
    :param seq: sekwencja do zakodowania
    """
    
    #modiffication to allow for Ns (by Magda)
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


def scan_sequences(sequences_file, filters_file):    
    """
    Funkcja skanujaca filtrem sekwencje. Tworzy splot wartosci w filtrze
    z macierza one-hot encoding reprezentujaca sekwencje. Zwraca slownik,
    w ktorym kluczem jest nazwa filtra, a wartoscia lista najlepszych podsekwencji,
    gdzie kazda z nich pochodzi od innej sekwencji (dla kazdej sekwencji wybor
    najlepszej podsekwencji)
    
    :param sequences_file: plik z sekwencjami ze zbioru treningowego
    :param filters_file: plik z wartosciami filtra
    """
    
    # slownik, w ktorym kluczem sa nazwy kolejnych filtrow,
    # a wartoscia lista krotek, gdzie na zerowej pozycji jest 
    # najlepsza podsekwencja dla kazdej sekwencji, a na pierwszej
    # wartosc splotu tej podsekwencji z filtrem
    best_sequences = {}
    
    sequences = read_sequences(sequences_file)
    
    network_name = os.path.basename(filters_file).replace("_filter.txt", "")
    df_above = pd.read_csv("./statistics/{}_filters_above_tresh.csv".format(network_name))
    filters_above_tresh = df_above["Filter"].tolist()
    
    nucleotides = {0: "A", 1: "C", 2: "G", 3: "T"}

    
    for sequence in sequences:
        
        # kodowanie sekwencji na one hot encoding
        matrix = one_hot_encode(sequence).transpose()
        
        with open(filters_file, "r") as file:
            lines = file.readlines()
            f = 0
            for i in range(1,len(lines)-1,5):
                # wczytanie filtra
                filter = np.loadtxt(lines[i:i+4])
                best = float('-inf')
                if "filter_{}".format(f) in filters_above_tresh:
                    # przechodzenie po macierzy danej sekwencji filtrem o rozmiarze (4,19)
                    for j in range(np.shape(matrix)[1]):
                        if np.shape(matrix[:,j:j+np.shape(filter)[1]])[1] == np.shape(filter)[1]:
                            # podmacierz sekwencji o rozmiarze (4,19)
                            submatrix = matrix[:,j:j+np.shape(filter)[1]]
                            # sprawdzenie ktora podsekwencja pomnozona przez wartosci
                            # filtra da najwieksza wartosc
                            convolution = filter*submatrix
                            total = np.sum(convolution)
                            if total > best:
                                best = total
                                best_seq = submatrix
                                where = np.argwhere(best_seq)
                                seq = "*" * np.shape(best_seq)[1]
                                # odtwarzanie sekwencji na podstawie macierzy one hot encoding
                                for index in where:
                                    seq = seq[:index[1]] + nucleotides[index[0]] + seq[index[1]+1:]
                    if "filter_"+str(f) in best_sequences.keys():
                        tmp = best_sequences["filter_"+str(f)]
                        tmp.append((seq, float(best)))
                    else:
                        best_sequences["filter_"+str(f)] = [(seq, float(best))]
                f += 1
    return best_sequences


def choose_best_subseq(dictionary, n=100):
    """
    Funkcja sluzaca do odnajdywania n najlepszych podsekwencji 
    w slowniku, w ktorym kluczem jest nazwa filtra, a wartoscia
    lista krotek, gdzie na zerowej pozycji jest najlepsza podsekwencja 
    dla kazdej sekwencji, ktora zostala przeskanowana tym filtrem,
    a na pierwszej pozycji wartosc splotu tej podsekwencji z filtrem.
    Ta funkcja sposrod tych wszystkich znalezionych podsekwencji wybierze 
    n najlepszych, ktore posluza do zbudowania macierzy PFM, a potem 
    PPM, a nastepnie PWM.
    Funkcja zwraca slownik, w ktorym kluczem jest nazwa filtra, a wartoscia
    lista najlepszych n podsekwencji
    
    :param dictionary: slownik, w ktorym kluczem jest nazwa filtra,
    a wartoscia lista najlepszych podsekwencji wraz z ich wartoscia splotu
    :param n: liczba najlepszych podsekwencji, ktore nalezy wybrac
    """

    dictionary_new = {}
    for key in dictionary:
        sorted_list = sorted(dictionary[key], key=lambda tup: tup[1], reverse=True)
        tmp = []
        for i in range(n):
            tmp.append(sorted_list[i][0])
        dictionary_new[key]=tmp

    return dictionary_new


def create_pfm_matrices(dictionary, dataset, network):
    """
    Funkcja tworzaca macierz PFM na podstawie slownika, w ktorym
    kluczem jest nazwa filtra, a wartoscia lista najlepszych
    n podsekwencji. Zapisuje utworzone macierze PFM dla kazdego
    filtra do pliku
    
    :param dictionary: slownik, w ktorym kluczem jest nazwa filtra, 
    a wartoscia lista najlepszych n podsekwencji
    """
    
    path = "./pfm_filters/{}/".format(dataset)
    if not os.path.exists(path):
        os.makedirs(path)
    
    with open(path+"{}_pfm.txt".format(network), "w") as fileout:
        fileout.write("MEME version 4\n\nALPHABET= ACGT\n\nstrands: + -\n\nBackground letter frequencies\nA 0.25 C 0.25 G 0.25 T 0.25\n\n")
        # dla kazdego filtra
        for filter in dictionary:
            fileout.write("MOTIF {}\nletter-probability matrix: alength= 4 w= 19\n".format(str(filter)))
            instances = []
            # tworzenie instancji Seq dla kazdej podsekwencji 
            for seq in dictionary[filter]:
                instances.append(Seq(seq))
            # tworzenie macierzy pfm dla zbioru podsekwencji dla kazdego filtra
            m = motifs.create(instances)
            pfm_matrix = m.counts
            # zamiana instancji Motif na macierz numpy
            # potrzebne, aby macierz pfm byla o wymiarach 19x4 a nie 4x19
            pfm_list = []
            for nucleotide in ["A", "C", "G", "T"]:
                pfm_list.append(list(pfm_matrix[nucleotide]))
            pfm_list = np.matrix(pfm_list).transpose()
            np.savetxt(fileout, pfm_list)
            fileout.write("\n")

    
    

def main():
    
    datasets = ["nonpromoter_active", "nonpromoter_inactive", "promoter_active", "promoter_inactive"]
    networks = ["alt-again-1", "alt2", "custom1", "patient_specific_thresh2_40000"]
    
    for dataset in datasets:
        print(dataset)
        for network in networks:
            print(network)
            scanned_seqs = scan_sequences("../Magisterka/dataset3/{}_10000.fa".format(dataset), "./filters/{}_filter.txt".format(network))
            #print(scanned_seqs, "\n")
            best_subseq = choose_best_subseq(scanned_seqs)
            #print(best_subseq, "\n")
            create_pfm_matrices(best_subseq, dataset, network)

if __name__ == "__main__":
    main()