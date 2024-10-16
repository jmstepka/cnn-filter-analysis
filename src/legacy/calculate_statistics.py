import os
import pandas as pd

def calculate_kmers_statistics(k, strict_intersection, networks_datasets_dict):
    all_kmers_dfs = {}
    kmers_perc_dict = {"Network" : [], "Dataset" : [], "% filters high IC" : [], "% kmers high IC" : [], 
                            "% filters TF" : [], "% kmers TF" : [], "% kmers high IC no TF" : []}
    for network, datasets in networks_datasets_dict.items():
        all_kmers_dfs[network] = {}
        for dataset in datasets:
            all_kmers_dict = {"Filter" : [], "Filter kmer" : []}
            
            #read all input filters and prepare rows in k-mers df
            input_directory = "./tomtom/{}/{}/tomtom_input/".format(dataset, network)
            for i_filename in os.listdir(input_directory):
                if "filter" in i_filename:
                    filter_name = i_filename.split(".")[0]
                    for i in range(19-k+1):
                        all_kmers_dict["Filter"].append(filter_name)
                        all_kmers_dict["Filter kmer"].append("{}_{}_{}".format(filter_name, i, i + k))
            all_kmers = pd.DataFrame(all_kmers_dict)
            
            all_kmers_in_analysis = pd.DataFrame()
            #read dfs with tomtom results (for all kmers, also low IC)
            directory = "./tomtom/{}/{}/tomtom_analysis/count_tf_hits/".format(dataset, network)
            for filename in os.listdir(directory):                 
                f = os.path.join(directory, filename)
                if os.path.isfile(f):
                    #dataframe with all kmers from onr filter
                    tf_hits = pd.read_csv(f)

                    #check if 7mer has low IC filter_112.csv
                    below_path = "./tomtom/{}/{}/kmers_below_tresh/{}".format(dataset, network, filename)
                    below_df = pd.read_csv(below_path)
                    tf_hits["High IC"] = ~tf_hits["Filter kmer"].isin(below_df["Filter"])
                    tf_hits["High IC"] = tf_hits["High IC"].astype(int)

                    #check if 7mer contains central nucleotide of DNAshape correlating 5mer

                    all_kmers_in_analysis = pd.concat([all_kmers_in_analysis, tf_hits]) 
                    
            all_kmers = pd.merge(all_kmers, all_kmers_in_analysis, on = "Filter kmer", how = "left").fillna(0)

            #replace counts of motifs above 1 with 1 
            all_kmers["Similarity to TF"] = all_kmers['Number of TF motifs found'].mask(all_kmers['Number of TF motifs found'] >= 1, 1)
            all_kmers_dfs[network][dataset] = all_kmers

            #count percentages
            kmers_perc_dict["Network"].append(network)
            kmers_perc_dict["Dataset"].append(dataset)
            filters_num = len(all_kmers["Filter"].unique())
            kmers_num = len(all_kmers["Filter kmer"].unique())
            #calculate sumr foe each column e.g. how many kmers have high IC
            filters_kmers_counts = all_kmers.groupby("Filter").sum()

            #"% filters high IC"
            perc_filters_high_ic = len(filters_kmers_counts[filters_kmers_counts["High IC"] > 0]) * 100 / filters_num
            kmers_perc_dict["% filters high IC"].append(perc_filters_high_ic)
            #"% 7mers high IC"
            perc_kmers_high_ic = len(all_kmers[all_kmers["High IC"] > 0]) * 100 / kmers_num
            kmers_perc_dict["% kmers high IC"].append(perc_kmers_high_ic)

            #"% filters TF"
            perc_filters_tf = len(filters_kmers_counts[filters_kmers_counts["Similarity to TF"] > 0]) * 100 / filters_num
            kmers_perc_dict["% filters TF"].append(perc_filters_tf)
            #"% 7mers TF"
            perc_7mers_tf = len(all_kmers[all_kmers["Similarity to TF"] > 0]) * 100 / kmers_num
            kmers_perc_dict["% kmers TF"].append(perc_7mers_tf)

            #"% 7mers high IC no TF"
            perc_highIC_no_TF = len(all_kmers[(all_kmers["Similarity to TF"] == 0) & (all_kmers["High IC"] == 1)]) * 100 / kmers_num
            kmers_perc_dict["% kmers high IC no TF"].append(perc_highIC_no_TF)


    kmers_perc_df = pd.DataFrame(kmers_perc_dict)
    return kmers_perc_df, all_kmers_dfs