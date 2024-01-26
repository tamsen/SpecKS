
#eventually this should be read in from an xml config or similar

class SpecKS_config:

    output_folder_root="/home/tamsen/Data/SpecKS_output"
    #output_folder_root="/Users/tamsen/Data/SpecKS_output"
    output_folder="TBD"
    dup_rate_parameters = (4, 2460)
    loss_rate_parameters = (4, 2053)
    num_replicates_per_gene_tree=5#10
    num_codons=30#1000
    tree_length=5
    path_to_sagephy = "/home/tamsen/Apps/sagephy/sagephy-1.0.0.jar"
    #path_to_sagephy = "/Users/tamsen/Apps/sagephy/sagephy-1.0.0.jar"
    # species_tree='(O:500,(P1:200,P2:200):300);'
    # dup_rate=0.00162
    # loss_rate=0.00194

    max_ks_for_hist_plot=2
    max_y_for_hist_plot=False
