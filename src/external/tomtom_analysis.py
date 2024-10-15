import os
import subprocess
import pandas as pd
from pathlib import Path
import csv
import shutil

def tomtom(kmers_files, dataset, network, hocomoco_path, thresh=0.5):
    """
    Funkcja tworzaca analize Tomtom. Wczytuje plik z kmerami
    z macierzy PPM i zapisuje analize Tomtom do folderow
    odpowiadajacym filtrom
    
    :param kmers_files: plik z kmerami z macierzy PPM
    :param dataset: analizowana klasa
    :param network: analizowana siec
    """

    path = "./tomtom/{}/{}/tomtom_output/".format(dataset, network) 

    if not os.path.exists(path):
        os.makedirs(path)
    else:
        shutil.rmtree(path)
        os.makedirs(path)
    
    for entry in os.scandir(kmers_files):
        if entry.path.endswith(".txt") and entry.is_file():
            filter = os.path.basename(entry.path).replace(".txt", "")
            subprocess.run("tomtom -incomplete-scores -dist kullback -o {} {} {} -thresh {} -verbosity 1".format(path+filter, entry.path, hocomoco_path, thresh), shell=True)

def count_found_motifs(dataset, network, k=7):
    """
    Funkcja tworzaca pliki csv, w ktorych znajduje sie analiza wynikow tomtom.
    Zlicza ile znaleziono motywow TF dla kazdego kmeru.
    Tworzy plik, w ktorym w wierszach znajduje sie kmer, a w kolumnach
    nazwy wszystkich TF. Na przecieciu stoi 1 gdy dany TF byl hitem dla tego
    kmeru, a 0 wpp
    
    :param dataset: analizowana klasa
    :param network: analizowana siec
    :param k: dlugosc kmeru
    """

    path = "./tomtom/{}/{}/tomtom_analysis/".format(dataset, network)
    
    if not os.path.exists(path):
        os.makedirs(path)
    
    path2 = "./tomtom/{}/{}/tomtom_analysis/count_tf_hits/".format(dataset, network)
    
    if not os.path.exists(path2):
        os.makedirs(path2)

    # lista wszystkich znalezionych TF dla wszystkich kmerow w kazdym filtrze
    tf_names = []
    
    # slownik, w ktorym kluczem jest nazwa kmeru, a wartoscia lista znalezionych dla niego 
    # motywow TF
    dictionary = {}
    # przechodzenie po katalogach odpowiadajacym filtrom
    for subdir, dirs, files in os.walk("./tomtom/{}/{}/tomtom_output/".format(dataset, network)):
        for file in files:
            if file.endswith(".tsv"):
                df = pd.read_csv(os.path.join(subdir, file), sep='\t')
                for index, row in df.iterrows():
                    if "H11" in str(row["Target_ID"]):
                        if row["Query_ID"] not in dictionary:
                            dictionary[row["Query_ID"]] = [row["Target_ID"]]
                        else:
                            dictionary[row["Query_ID"]].append(row["Target_ID"])
                            
                # tworzenie plikow zliczen znalezionych TF dla kazdego kmeru
                with open(path2+Path(subdir).stem+".csv", "w") as count_hits:
                    writer = csv.writer(count_hits)
                    writer.writerow(["Filter kmer", "Number of TF motifs found"])
                    for i in range(0, 19-k+1):
                        kmer = "{}_{}_{}".format(Path(subdir).stem, i, i+k)
                        if kmer in list(dictionary.keys()):
                            writer.writerow([kmer, len(dictionary[kmer])])
                        else:
                            writer.writerow([kmer, 0])
                            
                for key in dictionary:
                    for value in dictionary[key]:
                        if value not in tf_names:
                            tf_names.append(value)
    # tworzenie pliku przedstawiajacego ktore motywy TF wystapily 
    # w tym kmerze (w komorce 1 jesli spelnione, 0 wpp)       
    with open(path+ "tf_motifs_in_kmers.csv", "w") as fileout:
        writer = csv.writer(fileout)
        row = [i for i in tf_names]
        row.insert(0, "K-mer")
        writer.writerow(row)
        for key in dictionary:
            row2 = [key]
            for r in row:
                if r != "K-mer":
                    if r in dictionary[key]:
                        row2.append(1)
                    else:
                        row2.append(0)
            writer.writerow(row2)
                        

def aggregate_filters(dataset, network):
    """
    Funkcja tworzaca agregacje analizy Tomtom do filtrow.
    Pobiera hity dla pliku csv dla kmerow i tworzy nowy plik csv,
    w ktorym te dane sa agregowane do filtrow.
    
    :param dataset: analizowana klasa
    :param network: analizowana siec
    """
    
    path = "./tomtom/{}/{}/tomtom_analysis/".format(dataset, network)
    
    if not os.path.exists(path):
        os.makedirs(path)
    
    # wczytanie pliku, w ktorym w wierszach znajduje sie kmer, a w kolumnach nazwy wszystkich TF. 
    # Na przecieciu stoi 1 gdy dany TF byl hitem dla tego kmeru, a 0 wpp
    df = pd.read_csv("./tomtom/{}/{}/tomtom_analysis/tf_motifs_in_kmers.csv".format(dataset, network))
    
    # slownik, w ktorym kluczem jest nazwa filtra, a wartoscia lista motywow TF dla niego znalezionych
    dictionary = {}
    # nazwy TF z pliku z kmerami
    columns = list(df.keys())
    for i in range(300):
        # analizowany filtr
        filter = "filter_{}".format(i)
        # dla kazdego kmeru z pliku
        for index, row in df.iterrows():
            # wczytanie nazwy kmeru z pliku
            filter_tmp = row['K-mer'].split("_")
            # dla jakiego filtra byl to kmer
            filter_tmp_2 = "_".join(filter_tmp[:2])
            # jesli analizowany filtr byl w pliku z kmerami
            if filter == filter_tmp_2:
                # sprawdzenie jakie motywy TF znalezione dla kmeru tego filtra
                indices = [i for i, x in enumerate(row) if x == 1]
                # dodanie nazw TF znalezionych dla analizowanego filtra
                if filter in dictionary:
                    for ind in indices:
                        dictionary[filter].append(columns[ind])
                else:
                    dictionary[filter] = [columns[ind] for ind in indices]
    # unikanie powtorzen nazw TF dla filtra
    for key in dictionary:
        dictionary[key] = list(set(dictionary[key]))
    

    # zapisanie do pliku z agregacja do filtra
    with open(path+"tf_motifs_in_filters.csv", "w") as fileout:
        writer = csv.writer(fileout)
        row = [col for col in columns]
        row[0] = "Filter"
        row.insert(len(row), "Sum")
        writer.writerow(row)
        for key in dictionary:
            #values = [key]
            # dla kazdego filtra tworzenie 1 i 0 odpowiadajace wystapieniom motytow TF
            values = []
            # dla kazdego TF z pliku z kmerami
            for col in columns:
                if "K-mer" != col:
                    # jesli dany motyw TF byl znaleziony dla kmeru danego filtra
                    if col in dictionary[key]:
                        values.append(1)
                    # wpp
                    else:
                        values.append(0)
            suma = sum(values)
            values.append(suma)
            values.insert(0, key)
            writer.writerow(values)

def annotate_filters(dataset, network, hocomoco_annotation):
    """
    Funkcja tworzaca agregacje analizy Tomtom do rodzin.
    Pobiera hity dla pliku csv z agregacji do filtrow i tworzy 
    nowy plik csv, w ktorym te dane sa agregowane do rodzin.
    
    :param dataset: analizowana klasa
    :param network: analizowana siec
    """

    path = "./tomtom/{}/{}/tomtom_analysis/".format(dataset, network)
    
    if not os.path.exists(path):
        os.makedirs(path)
    
    # slownik, w ktorym kluczem jest nazwa TF a 
    # wartoscia nazwa rodziny, do ktorej nalezy
    dict_annotate = {}
    df = pd.read_csv(hocomoco_annotation, sep="\t")
    for index, row in df.iterrows():
        #dict_annotate[row["Model"]] = [row["TF family"],row["TF subfamily"]]
        # dla niektorych TF nie ma informacji o rodzinie
        if type(row["TF family"]) == str:
            dict_annotate[row["Model"]] = row["TF family"]
        else:
            dict_annotate[row["Model"]] = "No TF family info for this motif"
    
    vals = list(dict_annotate.values())
    for v in vals:
        if type(v) != str:
            print(v)
    
    # wczytanie pliku z agregacja do filtrow
    df = pd.read_csv(path+"tf_motifs_in_filters.csv")
    df = df.drop(columns=['Sum'])
    # slownik, w ktorym kluczem jest nazwa filtra a wartoscia lista
    # rodzin TF, ktore byly znalezione dla tego filtra
    dictionary = {}
    # nazwy TF obecne w pliku z agregacja do filtrow
    columns = list(df.keys())
    # lista nazw rodzin tych TF, ktore byly w pliku z agregacja dla filtrow
    columns_annotated = [dict_annotate[col] for col in columns if col != "Filter"]
    columns_annotated = list(set(columns_annotated))
    # dla kazdego filtra
    for index, row in df.iterrows():
        # ktore motywy TF znalezione dla tego filtra
        indices = [i for i, x in enumerate(row) if x == 1]
        # dodanie do slownika informacji o rodzinach TFs dla kazdego filtra
        if row["Filter"] in dictionary:
            for ind in indices:
                dictionary[row["Filter"]].append(dict_annotate[columns[ind]])
        else:
            dictionary[row["Filter"]] = [dict_annotate[columns[ind]] for ind in indices]
            
    # usuniecie powtorzen nazwy rodzin dla danego filtra
    for key in dictionary:
        dictionary[key] = list(set(dictionary[key]))
    # zapisanie do pliku z agregacja do rodzin
    with open(path+"tf_families_in_filters.csv", "w") as fileout:
        writer = csv.writer(fileout)
        row = [col for col in columns_annotated]
        row.insert(0, "Filter")
        row.insert(len(row), "Sum")
        writer.writerow(row)
        for key in dictionary:
            values = []
            for col in columns_annotated:
                if col != "Filter":
                    if col in dictionary[key]:
                        values.append(1)
                    else:
                        values.append(0)
            suma = sum(values)
            values.append(suma)
            values.insert(0, key)
            writer.writerow(values)


def main():
    
    datasets = ["nonpromoter_active", "nonpromoter_inactive", "promoter_active", "promoter_inactive"]
    networks = ["alt-again-1", "alt2", "custom1", "patient_specific_thresh2_40000"]
    
    #datasets = ["nonpromoter_active"]
    #networks = ["alt-again-1"]

    for dataset in datasets:
        print(dataset)
        for network in networks:
            print(network)
            # analiza tomtom
            tomtom("./tomtom/{}/{}/tomtom_input/".format(dataset, network), dataset, network)
            
            # zliczenia hitow i podsumowanie znalezionych TF dla kmerow
            count_found_motifs(dataset, network)
            aggregate_filters(dataset, network)
            annotate_filters(dataset, network)

if __name__ == "__main__":
    main()