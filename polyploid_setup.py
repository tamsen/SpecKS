import os


def make_polyploids(config):

    tree_1=setup_allopolyploid(config)
    tree_2=setup_autopolyploid(config)
    return [tree_1,tree_2]
def setup_allopolyploid(config):

    subfolder = config.output_folder
    SPC_time_MYA = 300
    WGD_time_MYA = 200
    species_name = "Allopolyploid_1"
    species_tree_result = polyploid_data(subfolder, species_name, SPC_time_MYA, WGD_time_MYA, config)
    return species_tree_result

def setup_autopolyploid(config):

    subfolder = config.output_folder
    SPC_time_MYA = 200
    WGD_time_MYA = 200
    species_name = "Autopolyploid_1"
    species_tree_result = polyploid_data(subfolder, species_name, SPC_time_MYA, WGD_time_MYA, config)
    return species_tree_result


class polyploid_data():

    simple_newick=""
    species_name=""
    SPC_time_MYA=0
    WGD_time_MYA=0
    FULL_time_MYA=0
    mode="TBD"
    general_sim_config=False
    species_subfolder=""
    subtree_subfolder=""
    analysis_step_num=1       #this will be

    # update at run-time
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