import numpy as np
import os
import csv

def create_pseudo_counts_matrix(pfm_matrix, epsilon = 1): 
    """
    Funkcja tworzaca macierz pseudo counts na podstawie macierzy pfm
    
    :param pfm: macierz pfm
    :param epsilon: wartosc pseudo counts
    """
    
    pseudo_counts_matrix = np.empty([pfm_matrix.shape[0], pfm_matrix.shape[1]])
    
    for (x,y), value in np.ndenumerate(pfm_matrix):
        pseudo_counts_matrix[x][y] = (value+(epsilon/4))/(100+epsilon)
    
    return pseudo_counts_matrix


def create_ppm_matrix(pseudo_counts_matrix):
    """
    Funkcja tworzaca macierz ppm na podstawie macierzy pfm z pseudo zliczeniami.

    :param pseudo_counts_matrix: macierz pfm z pseudo zliczeniami
    """
    
    ppm_matrix = np.empty([pseudo_counts_matrix.shape[0], pseudo_counts_matrix.shape[1]])
    for (x,y), value in np.ndenumerate(pseudo_counts_matrix):
        total = np.sum(pseudo_counts_matrix[x])
        ppm_matrix[x][y] = value/total
    return ppm_matrix
   

def create_pwm_matrix(pseudo_count_matrix):
    """
    Funkcja tworzaca macierz PWM na podstawie macierzy 
    pseudo counts (PPM with pseudocounts)
    
    :param pseudo_count_matrix: macierz PPM z pseudocounts
    """
    
    pwm_matrix = np.empty([pseudo_count_matrix.shape[0], pseudo_count_matrix.shape[1]])
   
    for (x,y), value in np.ndenumerate(pseudo_count_matrix):
        pwm_matrix[x][y] = np.log2(value/0.25) # 0.25 - czestosc tla
        
    total = np.sum(pwm_matrix, axis=1)
    for (x,y), value in np.ndenumerate(pwm_matrix):
        pwm_matrix[x][y] = (pwm_matrix[x][y])/total[y]

    return pwm_matrix
    

def ic_content(ppm_matrix):
    """
    Funkcja obliczaca ic content w macierzy ppm
    
    :param ppm_matrix: macierz ppm
    """
    
    ic_matrix = np.empty([ppm_matrix.shape[0], ppm_matrix.shape[1]])
    
    
    for r in range(len(ppm_matrix)):
        uncertainty = 0
        for value in ppm_matrix[r]:
            uncertainty += value * np.log2(value)
        uncertainty = (-1)*uncertainty
        ic_final = 2 - uncertainty
        row = []
        for value in ppm_matrix[r]:
            row.append(value*ic_final)
        ic_matrix[r] = row
    
    total = np.sum(ic_matrix)
    
    return total
    

def calculate_ic(pfm_file, dataset, network, k=7, ic_tresh = 6.5):
    """
    Funkcja sluzaca do obliczania IC dla macierzy PFM. 
    Z macierzy PFM tworzy pseudo counts a nastepnie PPM.
    Dzieli macierz PPM na kmery, dla ktorych oblicza IC.
    Jesli IC >= ic_tresh, zapisuje dany kmer (ppm) do pliku, 
    ktory bedzie plikiem wejsciowym do analizy tomtom.
    Jesli IC < ic_tresh, zapisuje wartosc IC dla tego kmera 
    do innego pliku (below_tresh)
    
    :param pfm_file: plik txt z macierza pfm reprezentujaca filtr
    :param dataset: rozpatrywana klasa
    :param network: rozpatrywana siec
    :param k: dlugosc kmeru
    :param ic_tresh: wartosc progu IC
    """
    
    path_input = "./tomtom/{}/{}/tomtom_input/".format(dataset, network)
    path_kmers = "./tomtom/{}/{}/kmers_below_tresh/".format(dataset, network)
    path_ic = "./statistics/ic_stats/{}_{}.csv".format(network, dataset) 
        
    if not os.path.exists(path_input):
        os.makedirs(path_input)

    if not os.path.exists(path_kmers):
        os.makedirs(path_kmers)
    
    with open(pfm_file) as file:
        lines_tmp = file.readlines()
        lines = lines_tmp[9:]
        for l in lines:
            if "letter-" in l:
                lines.remove(l)
            if l == "\n":
                lines.remove(l)
        ic_csv = open(path_ic, "w")
        writer_ic = csv.writer(ic_csv)
        writer_ic.writerow(["Filter", "IC"])
        for i in range(0,len(lines)-1,20):
            # wagi danego filtra
            weights = np.loadtxt(lines[i+1:i+20])
            # nazwa danego filtra
            filter = lines[i].split(" ")[1].strip()
            
            with open(path_input+filter+".txt", "w") as filterkmers:
                filterkmers.write("MEME version 4\n\nALPHABET= ACGT\n\nstrands: + -\n\nBackground letter frequencies\nA 0.25 C 0.25 G 0.25 T 0.25\n\n")
                with open(path_kmers+filter+".csv", "w") as belowtresh:
                    writer = csv.writer(belowtresh)
                    writer.writerow(["Filter", "IC"])
                    # szukanie kmerow
                    for i in range(0, 19-k+1):
                        # kmer wag filtra
                        kmer_freq = weights[i:(i+k)]
                        kmer_pseudo_counts = create_pseudo_counts_matrix(kmer_freq)
                        kmer_ppm = create_ppm_matrix(kmer_pseudo_counts)
                        info_content = ic_content(kmer_ppm)
                        if info_content >= ic_tresh:
                            #kmer_pwm = create_pwm_matrix(kmer_pseudo_counts)
                            filterkmers.write("MOTIF {}_{}_{}\nletter-probability matrix: alength= 4 w= {}\n".format(filter, i, i+k, k))
                            np.savetxt(filterkmers, kmer_ppm)
                            filterkmers.write("\n")
                            #print(kmer_ppm)
                        else:
                            writer.writerow(["{}_{}_{}".format(filter, i, i+k), info_content])
                        writer_ic.writerow(["{}_{}_{}".format(filter, i, i+k), info_content])
                            

def main():
    
    datasets = ["nonpromoter_active", "nonpromoter_inactive", "promoter_active", "promoter_inactive"]
    networks = ["alt-again-1", "alt2", "custom1", "patient_specific_thresh2_40000"]
    
    #datasets = ["nonpromoter_active"]
    #networks = ["alt-again-1"]

    for dataset in datasets:
        print(dataset)
        for network in networks:
            print(network)
            calculate_ic("./pfm_filters/{}/{}_pfm.txt".format(dataset, network), dataset, network)

if __name__ == "__main__":
    main()