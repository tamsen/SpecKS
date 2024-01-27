


def make_trees(config):

    tree_1=setup_allopolyploid(config.output_folder)
    tree_2=setup_autopolyploid(config.output_folder)
    return [tree_1,tree_2]
def setup_allopolyploid(subfolder):

    SPC_time_MYA = 300
    WGD_time_MYA = 250
    species_name = "Allopoly_1"
    species_tree_result = species_tree_data(subfolder, species_name, SPC_time_MYA, WGD_time_MYA)
    return species_tree_result

def setup_autopolyploid(subfolder):

    SPC_time_MYA = 300
    WGD_time_MYA = 300
    species_name = "Autopoly_1"
    species_tree_result = species_tree_data(subfolder, species_name, SPC_time_MYA, WGD_time_MYA)
    return species_tree_result


class species_tree_data():
    simple_newick=""
    species_name=""
    SPC_time=0
    WGD_time=0
    subfolder=0
    mode="TBD"
    def __init__(self, subfolder,species_name, SPC_time,WGD_time):

        self.species_tree_subfolder=os.path.join(subfolder,species_name)
        self.SPC_time = SPC_time
        self.WGD_time = WGD_time

        if (SPC_time==WGD_time):
            self.mode="AUTO"
        else:
            self.mode="AUTO"