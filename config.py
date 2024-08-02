# eventually this should be read in from an xml config or similar
import math
import xml.etree.ElementTree as ET


class SpecKS_config:

    output_folder_root = "/home/SpecKS_output"
    full_sim_time = 100 #MY

    mean_gene_birth_rate = 0.001359 #Ya-Long Guo reference
    mean_SSD_life_span= 5.77 #a half-life of 4 MY = mean life span of 5.77078016 ML
    mean_WGD_life_span= 44.72 #a half-life of 31 MY = mean life span of 44.723
    num_gene_trees_per_species_tree = 1000
    num_replicates_per_gene_tree = 0
    num_codons = 1000
    Ks_per_Myr = 0.01
    per_site_evolutionary_distance = 0.01268182
    evolver_random_seed = 137

    max_ks_for_hist_plot = 5
    max_y_for_hist_plot = False
    params_for_polyploids = []

    # debugging options
    stop_at_step = 999
    log_file_name = "log.txt"
    include_visualizations = False

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

            if (incoming_tag == "SpeciesTree"):
                for inner_layer in top_layer:
                    incoming_txt = inner_layer.text.strip()
                    incoming_tag = inner_layer.tag.strip()
                    if (incoming_tag == "full_sim_time"):
                        self.full_sim_time = int(incoming_txt)
                    if (incoming_tag == "polyploid"):
                        new_params = PolyploidParams(-1, -1, [],"foo")
                        for poly_layer in inner_layer:
                            incoming_txt = poly_layer.text.strip()
                            incoming_tag = poly_layer.tag.strip()
                            if (incoming_tag == "DIV_time_MYA"):
                                new_params.DIV_time_MYA = float(incoming_txt)
                            if (incoming_tag == "SPC_time_MYA"):
                                new_params.DIV_time_MYA = float(incoming_txt)
                            if (incoming_tag == "WGD_time_MYA"):
                                new_params.WGD_time_MYA = float(incoming_txt)
                            if (incoming_tag == "gene_div_time_distribution_parameters"):
                                new_params.divergence_distribution_parameters_list.append(
                                    parse_comma_separated_values(incoming_txt))

                            if (incoming_tag == "name"):
                                new_params.name = incoming_txt
                        self.params_for_polyploids.append(new_params)

            if (incoming_tag == "GeneTree"):
                for inner_layer in top_layer:
                    incoming_txt = inner_layer.text.strip()
                    incoming_tag = inner_layer.tag.strip()

                    if (incoming_tag == "mean_gene_birth_rate_GpMY"):
                        self.mean_gene_birth_rate = parse_float_or_false(incoming_txt)
                    if (incoming_tag == "mean_SSD_life_span_MY"):
                        self.mean_SSD_life_span = parse_float_or_false(incoming_txt)
                    if (incoming_tag == "mean_WGD_life_span_MY"):
                        self.mean_WGD_life_span = parse_float_or_false(incoming_txt)
                    if (incoming_tag == "SSD_half_life_MY"):
                        self.mean_SSD_life_span = half_life_to_mean_life(incoming_txt)
                    if (incoming_tag == "WGD_half_life_MY"):
                        self.mean_WGD_life_span = half_life_to_mean_life(incoming_txt)
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
                    if (incoming_tag == "per_site_evolutionary_distance"):
                        self.per_site_evolutionary_distance = float(incoming_txt)
                    if (incoming_tag == "evolver_random_seed"):
                        self.evolver_random_seed = int(incoming_txt)
def parse_tuple_string(tuple_string):
    if tuple_string.upper() == "FALSE":
        return False
    else:
        splat = tuple_string.replace("(", "").replace(")", "").split(",")
        data = [float(s) for s in splat]
        return data

def parse_float_or_false(input_string):
    if input_string.upper() == "FALSE":
        return False
    else:
        return float(input_string)

def half_life_to_mean_life(input_string):

    half_life=parse_float_or_false(input_string)
    if not half_life:
        return False
    else:
        ln2=math.log(2)
        return half_life/ln2
def parse_comma_separated_values(input_string):
    if input_string.upper() == "FALSE":
        return False
    else:
        cleaned_string=input_string.replace("(", "").replace(")", "")
        splat = cleaned_string.split(",")
        return splat


class PolyploidParams:
    DIV_time_MYA = 0
    WGD_time_MYA = 0
    divergence_distribution_parameters_list = []
    name = False

    def __init__(self, DIV_time_MYA, WGD_time_MYA, divergence_distribution_parameters_list, name):
        self.DIV_time_MYA = DIV_time_MYA
        self.WGD_time_MYA = WGD_time_MYA
        self.divergence_distribution_parameters_list = divergence_distribution_parameters_list
        self.name = name

    def to_xml(self):
        s1="\t<name>{0}</name>".format(self.name)
        s2="\t<DIV_time_MYA>{0}</DIV_time_MYA>".format(self.DIV_time_MYA)
        s3="\t<WGD_time_MYA>{0}</WGD_time_MYA>".format(self.WGD_time_MYA)
        return "\n".join([s1,s2,s3])
