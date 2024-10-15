import torch
from torch import nn
import math
import numpy as np
import os
import csv
import pandas as pd

from typing import Tuple

class CustomNetwork(torch.nn.Module):

    def __init__(self, seq_len=2000, num_channels=[300, 200, 200], kernel_widths=[19, 11, 7], pooling_widths=[3, 4, 4],
                 num_units=[2000, 4], dropout=0.5):
        super(CustomNetwork, self).__init__()
        paddings = [int((w-1)/2) for w in kernel_widths]
        self.seq_len = seq_len
        self.dropout = dropout
        self.params = {
            'input sequence length': seq_len,
            'convolutional layers': len(num_channels),
            'fully connected': len(num_units),
            'number of channels': num_channels,
            'kernels widths': kernel_widths,
            'pooling widths': pooling_widths,
            'units in fc': num_units,
            'dropout': dropout

        }

        conv_modules = []
        num_channels = [1] + num_channels
        for num, (input_channels, output_channels, kernel, padding, pooling) in \
                enumerate(zip(num_channels[:-1], num_channels[1:], kernel_widths, paddings, pooling_widths)):
            k = 4 if num == 0 else 1
            conv_modules += [
                nn.Conv2d(input_channels, output_channels, kernel_size=(k, kernel), padding=(0, padding)),
                nn.BatchNorm2d(output_channels),
                nn.ReLU(),
                nn.MaxPool2d(kernel_size=(1, pooling), ceil_mode=True)
            ]
            seq_len = math.ceil(seq_len / pooling)
        self.conv_layers = nn.Sequential(*conv_modules)

        fc_modules = []
        self.fc_input = 1 * seq_len * num_channels[-1]
        num_units = [self.fc_input] + num_units
        for input_units, output_units in zip(num_units[:-1], num_units[1:]):
            fc_modules += [
                nn.Linear(in_features=input_units, out_features=output_units),
                nn.ReLU(),
                nn.Dropout(p=self.dropout)
            ]
        self.fc_layers = nn.Sequential(*fc_modules)

    def forward(self, x):
        x = self.conv_layers(x)
        x = x.view(-1, self.fc_input) 
        x = self.fc_layers(x)
        return torch.sigmoid(x)


def get_filters(modelfile: str, fileout: str, filter_shape: Tuple[int, int]):
    """
    Funkcja sluzaca do wyodrebnienia filtrow z pierwszej (zerowej)
    warstwy sieci i zapisujaca je do pliku
    
    :param modelfile: plik z wyuczonym modelem (siecia)
    :param fileout: plik do ktorego zapisane bada wartosci filtrow
    """
    
    path = './filters/'
    if not os.path.exists(path):
        os.makedirs(path)
    
    use_cuda = torch.cuda.is_available()
    device = torch.device("cuda:0" if use_cuda else "cpu")
    model = CustomNetwork()
    # wczytanie wytrenowanego modelu
    model.load_state_dict(torch.load(modelfile, map_location=torch.device(device)))
    for name, param in model.named_parameters():
        # zerowa warstwa konwolucyjna
        if name == "conv_layers.0.weight":
            with open(path+fileout, "w") as f:
                for i in range(len(param)):
                    f.write(">filter_"+str(i)+"\n")
                    filter = param[i].detach().numpy().reshape(filter_shape)
                    np.savetxt(f,filter)


def calculate_statistics_filters():
    """
    Funkcja obliczajaca statystyki dotyczace filtrow. Tworzy plik csv, 
    w ktorym znajduje sie informacja o max, min, mean, std i medianie filtra
    """
    
    path = "./statistics/"
    if not os.path.exists(path):
        os.makedirs(path)
        
    for entry in os.scandir("./filters/"):
        if entry.path.endswith(".txt") and entry.is_file():
            with open(entry.path, "r") as file:                 
                lines = file.readlines()
                f = 0
                # pobranie nazwy sieci i utworzenie nazwy pliku
                filename = path+os.path.basename(entry.path).replace("_filter.txt", "_stats.csv")
                with open(filename, 'w') as fileout:
                    writer = csv.writer(fileout)
                    writer.writerow(["Filter", "Max", "Min", "Mean", "Std", "Median"])
                    for i in range(1,len(lines)-1,5):
                        # wczytanie filtra
                        filter = np.loadtxt(lines[i:i+4])
                        max = np.max(filter)
                        min = np.min(filter)
                        mean = np.mean(filter)
                        std = np.std(filter)
                        median = np.median(filter)
                        fil = "filter_{}".format(f)
                        writer.writerow([fil, max, min, mean, std, median])
                        f += 1

def divide_filters_weights(filter_stats, thresh = 1e-5):
    """
    Funkcja sluzaca do podzielenia zbioru filtrow na te, ktore
    maja srednie wartosci wag powyzej progu i te, ktore maja 
    te wartosci ponizej. Zapisane sa one do oddzielnych plikow,
    odpowiednio: {siec}_filter_above_tresh.csv oraz
    {siec}_filter_below_tresh.csv
    
    :param filter_stats: plik z statystykami dla filtrow w danej sieci
    :param tresh: prog odciecia sredniej wartosci filtrow
    """

    path = "./statistics/"
    if not os.path.exists(path):
        os.makedirs(path)
    
    stats = pd.read_csv(filter_stats)
    
    below_tresh = []
    above_tresh = []
    
    for index, row in stats.iterrows():
        # jesli srednia wartosci w filtrze powyzej progu
        if abs(row["Mean"]) > thresh:
            above_tresh.append([row["Filter"], row["Mean"]])
        # wpp
        else:
            below_tresh.append([row["Filter"], row["Mean"]])
    
    df_below = pd.DataFrame(below_tresh, columns=["Filter", "Mean"])
    df_above = pd.DataFrame(above_tresh, columns=["Filter", "Mean"])
    
    filter = os.path.basename(filter_stats).replace("_stats.csv", "")
    
    df_below.to_csv(path+filter+"_filters_below_tresh.csv")
    df_above.to_csv(path+filter+"_filters_above_tresh.csv")


def main():
    
    networks = ["alt-again-1", "alt2", "custom1", "patient_specific_thresh2_40000"]
    
    '''
    # get filters from all networks
    for network in networks:
        get_filters("../Magisterka/results_for_magda/{}_last.model".format(network), "{}_filter.txt".format(network))
    
    # calculate filter statistics
    calculate_statistics_filters()
    '''
    
    # get filters the weights of which meet threshold
    for network in networks:
        divide_filters_weights("./statistics/{}_stats.csv".format(network))

if __name__ == "__main__":
    main()