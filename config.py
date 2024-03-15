# eventually this should be read in from an xml config or similar
import xml.etree.ElementTree as ET


class SpecKS_config:
    output_folder_root = "/home/tamsen/Data/SpecKS_output"
    # output_folder_root="/Users/tamsen/Data/SpecKS_output"
    path_to_sagephy = "/home/tamsen/Apps/sagephy/sagephy-1.0.0.jar"
    # path_to_sagephy = "/Users/tamsen/Apps/sagephy/sagephy-1.0.0.jar"

    full_sim_time = 500
    dup_rate_parameters = (4, 2460)
    loss_rate_parameters = (4, 2053)
    branch_relaxation_parameters = ["ACRY07", "1", "0.000001"]
    num_gene_trees_per_species_tree = 4  # 10
    num_replicates_per_gene_tree = 3  # 10
    num_codons = 10  # 1000
    Ks_per_Myr = 0.01

    max_ks_for_hist_plot = 5
    max_y_for_hist_plot = False
    params_for_polyploids = []

    # debugging options
    stop_at_step = 999
    log_file_name = "log.txt"
    include_visualizations = True

    def __init__(self, config_file):

        mytree = ET.parse(config_file)
        myroot = mytree.getroot()

        for top_layer in myroot:

            incoming_tag = top_layer.tag.strip()
            incoming_txt = top_layer.text.strip()
            if (incoming_tag == "StopAtStep"):
                self.stop_at_step = int(incoming_txt)
            if (incoming_tag == "LogFileName"):
                if incoming_txt.upper() == "FALSE":
                    self.log_file_name = False
                else:
                    self.log_file_name = incoming_txt
            if (incoming_tag == "IncludeVisuals"):
                if incoming_txt.upper() == "FALSE":
                    self.include_visualizations = False
                else:
                    self.include_visualizations = True

            if (incoming_tag == "Histogram"):
                for inner_layer in top_layer:
                    incoming_txt = inner_layer.text.strip()
                    incoming_tag = inner_layer.tag.strip()
                    if (incoming_tag == "maxKs"):
                        self.max_ks_for_hist_plot = int(incoming_txt)
                    elif (incoming_tag == "maxY"):
                        if incoming_txt.upper() == "FALSE":
                            self.max_y_for_hist_plot = False
                        else:
                            self.max_ks_for_hist_plot = int(incoming_txt)

            if (incoming_tag == "Paths"):
                for inner_layer in top_layer:
                    incoming_txt = inner_layer.text.strip()
                    incoming_tag = inner_layer.tag.strip()
                    if (incoming_tag == "output_folder_root"):
                        self.output_folder_root = incoming_txt
                    elif (incoming_tag == "path_to_sagephy"):
                        self.path_to_sagephy = incoming_txt

            if (incoming_tag == "SpeciesTree"):
                for inner_layer in top_layer:
                    incoming_txt = inner_layer.text.strip()
                    incoming_tag = inner_layer.tag.strip()
                    if (incoming_tag == "full_sim_time"):
                        self.full_sim_time = int(incoming_txt)
                    if (incoming_tag == "polyploid"):
                        new_params = PolyploidParams(-1, -1, "foo")
                        for poly_layer in inner_layer:
                            incoming_txt = poly_layer.text.strip()
                            incoming_tag = poly_layer.tag.strip()
                            if (incoming_tag == "SPC_time_MYA"):
                                new_params.SPC_time_MYA = int(incoming_txt)
                            if (incoming_tag == "WGD_time_MYA"):
                                new_params.WGD_time_MYA = int(incoming_txt)
                            if (incoming_tag == "name"):
                                new_params.name = incoming_txt
                        self.params_for_polyploids.append(new_params)

            if (incoming_tag == "GeneTree"):
                for inner_layer in top_layer:
                    incoming_txt = inner_layer.text.strip()
                    incoming_tag = inner_layer.tag.strip()
                    if (incoming_tag == "dup_rate_parameters"):
                        self.dup_rate_parameters = parse_tuple_string(incoming_txt)
                    if (incoming_tag == "loss_rate_parameters"):
                        self.loss_rate_parameters = parse_tuple_string(incoming_txt)
                    if (incoming_tag == "branch_relaxer_parameters"):
                        self.branch_relaxation_parameters = parse_branch_relaxation_string(incoming_txt)
                    if (incoming_tag == "num_gene_trees_per_species_tree"):
                        self.num_gene_trees_per_species_tree = int(incoming_txt)

            if (incoming_tag == "SequenceEvolution"):
                for inner_layer in top_layer:
                    incoming_txt = inner_layer.text.strip()
                    incoming_tag = inner_layer.tag.strip()
                    if (incoming_tag == "num_replicates_per_gene_tree"):
                        self.num_replicates_per_gene_tree = int(incoming_txt)
                    if (incoming_tag == "num_codons"):
                        self.num_codons = int(incoming_txt)
                    if (incoming_tag == "Ks_per_Myr"):
                        self.Ks_per_Myr = float(incoming_txt)


def parse_tuple_string(tuple_string):
    if tuple_string.upper() == "FALSE":
        return False
    else:
        splat = tuple_string.replace("(", "").replace(")", "").split(",")
        data = [float(s) for s in splat]
        return data


def parse_branch_relaxation_string(branch_relaxation_string):
    if branch_relaxation_string.upper() == "FALSE":
        return False
    else:
        splat = branch_relaxation_string.split(",")
        return splat


class PolyploidParams:
    SPC_time_MYA = 0
    WGD_time_MYA = 0
    name = False

    def __init__(self, SPC_time_MYA, WGD_time_MYA, name):
        self.SPC_time_MYA = SPC_time_MYA
        self.WGD_time_MYA = WGD_time_MYA
        self.name = name
