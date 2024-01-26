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

    step_num=1
    #Time since WGD: 5,10, 15,50,100,200 MYA. Total tree length 500 MY. Make allo and autopoly examples.
    print("{}. Make species trees (SaGePhy or by hand)".format(step_num))
    species_tree = species_tree_maker.make_species_trees(conf, step_num)

    step_num=2
    print("{}. Make gene trees (SaGePhy)".format(step_num))
    num_gene_trees=2#10
    gene_tree_results_by_tree_name = gene_tree_maker.run_sagephy(conf,
                                                                 species_tree, num_gene_trees,step_num)
    step_num = 3
    print("{0}. Relax gene trees (SaGePhy)".format(step_num))
    num_gene_trees=2#10
    relaxed_gene_tree_results_by_tree_name = gene_tree_relaxer.relax(conf,
                                                                 species_tree, num_gene_trees,step_num)
    step_num = 4
    print("{0}. Evolve sequences through gene trees (Evolver)".format(step_num))
    evolver_results_by_gene_tree=gene_evolver.run_evolver(conf, gene_tree_results_by_tree_name,step_num)

    step_num = 5
    print("{0}. Prune trees. at every time step, cull a certain percent of what remains")
    tree_pruner.prune_gene_trees(conf,step_num)

    step_num = 6
    print("{0}. Get Ks for trees (Codeml)".format(step_num))
    #Make 1000 gene trees per species tree. Tiley adds that "we simulated variation in the size of the gene families by altering the number of genes at the root of each gene tree by allowing the root age to be exponentially distributed with a mean of 600 Ma"
    codeml_results_by_replicate_num = ks_calculator.run_codeml(conf,evolver_results_by_gene_tree,step_num)

    step_num = 7
    print("{0}. Plot histograms (matplotlib, or Ks rates)".format(step_num))
    ks_histogramer.run_Ks_histogramer(conf, codeml_results_by_replicate_num,
                                      gene_tree_results_by_tree_name,step_num)


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
