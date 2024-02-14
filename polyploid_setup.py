import os


def make_polyploids(config):

    subfolder = config.output_folder
    polyploids=[]
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

        polyploid = polyploid_data(subfolder, species_name, SPC_time_MYA, WGD_time_MYA, config)
        polyploids.append(polyploid)

    return polyploids

class polyploid_data():

    simple_newick=""
    species_name=""
    SPC_time_MYA=0
    WGD_time_MYA=0
    FULL_time_MYA=0
    subgenome_names = ["P1", "P2"]
    original_genome_name = ["P"]
    mode="TBD"
    simulation_legs=[]
    general_sim_config=False
    species_subfolder=""      #this will be updated at run-time
    subtree_subfolder=""      #this will be updated at run-time
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
            first_simulation_leg=sim_time_interval_forward_in_time(0,
                                                                    self.FULL_time_MYA-self.WGD_time_MYA,
                                                                    self.original_genome_name)
            second_simulation_leg=sim_time_interval_forward_in_time(self.FULL_time_MYA-self.WGD_time_MYA,
                                                                    self.FULL_time_MYA,self.subgenome_names)
            self.simulation_legs=[first_simulation_leg,second_simulation_leg]
        else:
            self.mode="ALLO"
            single_simulation_leg=sim_time_interval_forward_in_time(0,
                                                                    self.FULL_time_MYA,self.subgenome_names)
            self.simulation_legs = [single_simulation_leg]

    def is_allo(self):

        if self.mode=="ALLO":
            return True
        else:
            return False

    def WGD_time_as_ks(self):

        if self.general_sim_config:
            return self.general_sim_config.Ks_per_Myr * self.WGD_time_MYA
        else:
            return False

    def SPEC_time_as_ks(self):

        if self.general_sim_config:
            return self.general_sim_config.Ks_per_Myr * self.SPC_time_MYA
        else:
            return False


class sim_time_interval_forward_in_time():

    interval_start_time_MY=0
    interval_end_time_MY=0
    subgenomes_during_this_interval =[]
    def __init__(self, start_time_MY, end_time_MY,subgenomes):
        self.interval_start_time_MY = start_time_MY
        self.interval_end_time_MY = end_time_MY
        self.subgenomes_during_this_interval = subgenomes

    def leg_length(self):
        return self.interval_end_time_MY-self.interval_start_time_MY