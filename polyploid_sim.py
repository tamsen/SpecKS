"""polyploid_sim.py
Simulates the genome evolution of a single polyploid, as a step-by-step process.
"""

import log
from pipeline_modules import ks_histogramer, ks_calculator, gene_evolver, \
    species_tree_maker, results_organizer, gene_tree_maker, custom_GBD_model, gene_shedder

from Bio import Phylo
from io import StringIO
def run_sim(polyploid):


    #the allopolyploid has only one simulation leg for the full sim
    only_simulation_leg=polyploid.simulation_legs[0]
    evolver_random_seed=polyploid.general_sim_config.evolver_random_seed
    polyploid_genomes_of_interest = ['P1', 'P2']
    genomes_of_interest_by_species={"outgroup":['O'],polyploid.species_name:polyploid_genomes_of_interest}
    all_genomes_of_interest= polyploid_genomes_of_interest + ['O']

    log.write_to_log("\n\n{0}. Make species trees (custom code)".format(polyploid.analysis_step_num))
    species_trees = species_tree_maker.make_species_trees(polyploid)
    if polyploid.analysis_step_num > polyploid.general_sim_config.stop_at_step:
        return

    log.write_to_log("\n\n{0}. Make gene trees given gradual speciation".format(polyploid.analysis_step_num))
    base_gene_tree_newicks_by_tree_name = gene_tree_maker.make_randomized_gene_trees(polyploid,species_trees[0])
    if polyploid.analysis_step_num > polyploid.general_sim_config.stop_at_step:
        return

    log.write_to_log("\n\n{0}. Custom GBD model".format(polyploid.analysis_step_num))
    gene_data_by_gt_name = custom_GBD_model.run_custom_GBD_model(polyploid,all_genomes_of_interest,
                                                                 only_simulation_leg,
                                                                 base_gene_tree_newicks_by_tree_name)
    if polyploid.analysis_step_num > polyploid.general_sim_config.stop_at_step:
        return

    log.write_to_log("\n\n{0}. Prune WGD genes that will be dead before the end of the sim".format(polyploid.analysis_step_num))
    pruned_gene_data_by_gt_name = gene_shedder.shed_genes(polyploid,gene_data_by_gt_name,"P1")
    if polyploid.analysis_step_num > polyploid.general_sim_config.stop_at_step:
        return

    log.write_to_log("\n\n{0}. Evolve sequences through gene trees (Evolver)".format(polyploid.analysis_step_num))
    evolver_results_by_gene_tree = gene_evolver.run_evolver(polyploid, pruned_gene_data_by_gt_name,
                                                            evolver_random_seed)
    if polyploid.analysis_step_num > polyploid.general_sim_config.stop_at_step:
        return

    log.write_to_log("\n\n{0}. Get Ks for trees (Codeml)".format(polyploid.analysis_step_num))
    Ks_results_by_species_by_replicate_num = ks_calculator.run_codeml(polyploid,
                                                               genomes_of_interest_by_species,
                                                               pruned_gene_data_by_gt_name,
                                                               evolver_results_by_gene_tree)
    if polyploid.analysis_step_num > polyploid.general_sim_config.stop_at_step:
        return

    log.write_to_log("\n\n{0}. Plot histograms (matplotlib)".format(polyploid.analysis_step_num))
    final_result_files_by_species_by_replicate = ks_histogramer.run_Ks_histogramer(polyploid,
                                                                        genomes_of_interest_by_species,
                                                                        Ks_results_by_species_by_replicate_num)
    if polyploid.analysis_step_num > polyploid.general_sim_config.stop_at_step:
        return

    log.write_to_log("\n\n{0}. Collate results".format(polyploid.analysis_step_num))
    results_organizer.collate_results(polyploid, final_result_files_by_species_by_replicate)

    log.write_to_log("\n\n" + polyploid.species_name + " complete.\n\n")


#utility for debugging code..not actually part of the execution path
def check_for_gt_disparity(gene_tree_results_by_tree_name):
    for gt_name, gene_tree_result in gene_tree_results_by_tree_name.items():
        test_tree = Phylo.read(StringIO(gene_tree_result.simple_newick), "newick")
        X1 = Phylo.to_networkx(gene_tree_result.tree)
        X2 = Phylo.to_networkx(test_tree)
        len_nodes1 = len(list(X1.nodes))
        len_nodes2 = len(list(X2.nodes))
        print("nodes1=\t" + str(len_nodes1))
        print("nodes2=\t" + str(len_nodes2))
        if len_nodes1 != len_nodes2:
            print("big problem here!!")
            raise ValueError("caught the bug")