import os


def make_polyploids(config):

    subfolder = config.output_folder
    species_tree_results=[]
    num_allos=0
    num_autos=0

    for params in  config.params_for_polyploids:
        SPC_time_MYA=params.SPC_time_MYA
        WGD_time_MYA=params.WGD_time_MYA

        if (SPC_time_MYA==WGD_time_MYA):
            num_autos=num_autos+1
            species_name="Autopolyploid_" + str(num_autos)
        else:
            num_allos=num_allos+1
            species_name="Allopolyploid_" + str(num_allos)

        species_tree_result = polyploid_data(subfolder, species_name, SPC_time_MYA, WGD_time_MYA, config)
        species_tree_results.append(species_tree_result)

    return species_tree_results

class polyploid_data():

    simple_newick=""
    species_name=""
    SPC_time_MYA=0
    WGD_time_MYA=0
    FULL_time_MYA=0
    subgenome_names = ["P1", "P2"]
    original_genome_name = ["P"]
    mode="TBD"
    general_sim_config=False
    species_subfolder=""
    subtree_subfolder=""
    analysis_step_num=1       #this will be updated at run-time
    def __init__(self, subfolder,species_name, SPC_time,WGD_time,general_sim_config):

        self.species_subfolder=os.path.join(subfolder, species_name)
        self.general_sim_config =general_sim_config
        self.SPC_time_MYA = SPC_time
        self.WGD_time_MYA = WGD_time
        self.FULL_time_MYA = general_sim_config.full_sim_time
        self.species_name = species_name

        if (SPC_time==WGD_time):
            self.mode ="AUTO"
        else:
            self.mode="ALLO"

    def is_allo(self):

        if self.mode=="ALLO":
            return True
        else:
            return False