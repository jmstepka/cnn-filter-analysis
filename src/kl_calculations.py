import numpy as np
import pandas as pd
import os

def add_pseudocount(freq, sites, pseudocount = 0.1):
    '''
    Add pseudocount assuming uniform background.
    '''
    return ((freq * sites) + (pseudocount * 0.25)) / (sites + pseudocount)

def calculate_kl(query_freq, target_freq):
    '''
    Calculate Kullback-Leibler divergence between two positions.
    '''
    avg_kull = ((query_freq * np.log10(query_freq/target_freq)) + (target_freq * np.log10(target_freq/query_freq))) / 2
    return avg_kull

def calculate_kl_distance():
    hocomoco_models_path = "data/HOCOMOCOv11_core_HUMAN_mono_meme_format.meme"
    alphabet_size = 4
    pseudocount = 0.1

    for species in ["hg38", "pt06", "rm10"]:
        print(species)
        for network_type in ["_5k_last", "_5k_again_last"]:  
            network = species + network_type
            print(network)
            for dataset in ["pa", "na", "pi", "ni"]:   
                print(dataset)
                results_dict = {"Species" : [],
                                "Network" : [],
                                "Dataset" : [],
                                "Query" : [],
                                "Target" : [],
                                "Offset" : [],
                                "KL distance" : []}

                tomtom_input_path =  "tomtom/{}_5k_{}/{}/tomtom_input/".format(species, dataset, network)

                for pfm_file in os.listdir(tomtom_input_path):
                    if ".txt" in pfm_file:
                        #print(pfm_file)
                        with open(os.path.join(tomtom_input_path,pfm_file)) as file:
                            motifs = file.read().split("\nMOTIF ")[1:]
                            for motif in motifs:
                                qname = motif.split("\n")[0]
                                matrix = np.loadtxt(motif.split("\n")[2:])
                                length = len(matrix)                 

                                #iterate over target motifs
                                with open(hocomoco_models_path) as tfile:
                                    tmotifs = tfile.read().split("\nMOTIF ")[1:]
                                    for tmotif in tmotifs:
                                        tname = tmotif.split("\n")[0]
                                        tmatrix = np.loadtxt(tmotif.split("\nURL")[0].split("\n")[2:])
                                        nsites = int(tmotif.split("\n")[1].split("nsites= ")[-1])
                                        tlength = len(tmatrix)
                                        
                                        #reverse complement target
                                        rc_tmatrix = np.flip(tmatrix)
                                            
                                        distances_dict = {}

                                        #sliding window on target motif  
                                        for iq in range(0, length):
                                            #we want to include only full alignments
                                            for it in range(iq, iq + tlength-length + 1):
                                                offset = it - iq                        
                                                distance = 0
                                                rc_distance = 0

                                                for ia in range(0, alphabet_size):                                                
                                                    #get frequencies for the analyzed position
                                                    q_freq = matrix[iq,ia]
                                                    t_freq = tmatrix[it,ia] 
                                                    rc_t_freq = rc_tmatrix[it,ia]
                                                    
                                                    #add pseudocount to frequencies in the target motif
                                                    #(motifs from hocomoco may contain zeroes)                                                
                                                    t_freq = add_pseudocount(t_freq, nsites) 
                                                    rc_t_freq = add_pseudocount(rc_t_freq, nsites)  
                                                    
                                                    #calculate KL distance
                                                    avg_kull = calculate_kl(q_freq, t_freq)
                                                    avg_kull_rc = calculate_kl(q_freq, rc_t_freq)
                                                    distance += avg_kull
                                                    rc_distance += avg_kull_rc
                                                
                                                #change KL dist to negative value for better visualization
                                                #with log scale
                                                distance = -1 * distance
                                                rc_distance = -1 * rc_distance
                                                
                                                #sum position scores to calculate KL dist for the whole alignment
                                                if offset in distances_dict.keys():
                                                    distances_dict[offset] += distance
                                                    distances_dict[-(offset+1)] += rc_distance
                                                else:
                                                    distances_dict[offset] = distance
                                                    distances_dict[-(offset+1)] = rc_distance  
                                                    
                                                #select the lowest value (take all in case of ties)
                                                offset_dist_list = [[o, d] for o,d in distances_dict.items()]
                                                offset_dist_list_sorted = sorted(offset_dist_list, key=lambda x: x[1], reverse=True)  
                                                lowest_score = offset_dist_list_sorted[0][1]
                                                lowest_scores_dict = {offset_dist_list_sorted[0][0] : offset_dist_list_sorted[0][1]}
                                                for o, d in offset_dist_list_sorted:
                                                    if d == lowest_score:
                                                        lowest_scores_dict[o] = d
                                                    else:
                                                        break                                                   
                                                    

                                        for offset in lowest_scores_dict.keys():
                                            results_dict["Species"].append(species)
                                            results_dict["Network"].append(network)
                                            results_dict["Dataset"].append(dataset)
                                            results_dict["Query"].append(qname)
                                            results_dict["Target"].append(tname)
                                            results_dict["Offset"].append(offset)
                                            results_dict["KL distance"].append(lowest_scores_dict[offset])
                results = pd.DataFrame(results_dict)
                results.to_csv("output_50/{}_{}_KLdist.csv".format(network, dataset), sep = "\t", index = False) 