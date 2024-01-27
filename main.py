import os

import species_tree_maker
import gene_tree_maker
import gene_tree_relaxer
import gene_evolver
import tree_pruner
import ks_calculator
import ks_histogramer
import config
from datetime import datetime


def run_sim():

    conf = setup()

    #Time since WGD: 5,10, 15,50,100,200 MYA. Total tree length 500 MY. Make allo and autopoly examples.
    print("\n\n{0}. Make species trees (custom code)".format(conf.sim_step_num))
    species_tree = species_tree_maker.make_species_trees(conf)

    #for species_tree  in species_trees:

    print("\n\n{0}. Make gene trees (SaGePhy)".format(conf.sim_step_num))
    gene_tree_results_by_tree_name = gene_tree_maker.run_sagephy(conf,species_tree)

    print("\n\n{0}. Relax gene trees (SaGePhy)".format(conf.sim_step_num))
    relaxed_gene_tree_results = gene_tree_relaxer.relax(conf,gene_tree_results_by_tree_name)

    print("\n\n{0}. Evolve sequences through gene trees (Evolver)".format(conf.sim_step_num))
    evolver_results_by_gene_tree=gene_evolver.run_evolver(conf, relaxed_gene_tree_results)

    print("\n\n{0}. Prune trees. ".format(conf.sim_step_num) +
            "At every time step post WGD, cull a certain percent of what remains. (custom code)")
    tree_pruner.prune_gene_trees(conf)

    print("\n\n{0}. Get Ks for trees (Codeml)".format(conf.sim_step_num))
    codeml_results_by_replicate_num = ks_calculator.run_codeml(conf,
                                    relaxed_gene_tree_results,evolver_results_by_gene_tree)

    print("\n\n{0}. Plot histograms (matplotlib)".format(conf.sim_step_num))
    ks_histogramer.run_Ks_histogramer(conf, codeml_results_by_replicate_num,relaxed_gene_tree_results)

    print("\n\nSpecKS complete")

def setup():

    now = datetime.now()
    date_time = now.strftime("m%md%dy%Y_h%Hm%Ms%S")
    conf = config.SpecKS_config()
    conf.output_folder = conf.output_folder_root + "_" + date_time

    print(conf.output_folder)
    if not os.path.exists(conf.output_folder):
        os.makedirs(conf.output_folder)

    print("Current environment:")
    print(str(os.environ))

    cwd=os.getcwd()
    print("Current Working Directory:\t" + cwd)

    return conf


if __name__ == '__main__':
    run_sim()

# See PyCharm help at https://www.jetbrains.com/help/pycharm/
