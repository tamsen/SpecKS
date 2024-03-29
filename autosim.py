import log
from pipeline_modules import ks_histogramer, ks_calculator, gene_evolver, \
    species_tree_maker, root_seq_maker, results_organizer, gene_tree_maker, custom_GBD_model, gene_shedder


def run_autosim(polyploid):


    #the autopolyploid has two simulation legs, one before speciation/wgd, and one after,
    preWGD_simulation_leg=polyploid.simulation_legs[0]
    postWGD_simulation_leg=polyploid.simulation_legs[1]
    polyploid_subgenomes = ['P1', 'P2']

    log.write_to_log("\n\n{0}. Make species trees (custom code)".format(polyploid.analysis_step_num))
    species_tree = species_tree_maker.make_species_trees(polyploid)
    if polyploid.analysis_step_num > polyploid.general_sim_config.stop_at_step:
        return

    log.write_to_log("\n\n{0}. Make gene trees given abrupt speciation".format(polyploid.analysis_step_num))
    base_gene_tree_newicks_by_tree_name = gene_tree_maker.make_randomized_gene_trees(polyploid,species_tree[0])
    if polyploid.analysis_step_num > polyploid.general_sim_config.stop_at_step:
        return

    log.write_to_log("\n\n{0}. Custom GBD model".format(polyploid.analysis_step_num))
    gene_data_by_gt_name = custom_GBD_model.run_custom_GBD_model(polyploid,preWGD_simulation_leg.subgenomes_during_this_interval,
                                                                 preWGD_simulation_leg,
                                                                 base_gene_tree_newicks_by_tree_name)
    if polyploid.analysis_step_num > polyploid.general_sim_config.stop_at_step:
        return

    # For the first leg of the autopolyploid sim, we evolve sequences from
    # what ever the start time was (say, 500 MYA) to the time of WGD (say, 300 MYA)

    first_leg_random_seed=137
    log.write_to_log("\n\n{0}. Evolve sequences through gene trees (Evolver)".format(polyploid.analysis_step_num))
    evolver_results_by_gene_tree = gene_evolver.run_evolver(polyploid, gene_data_by_gt_name,
                                                            first_leg_random_seed)
    if polyploid.analysis_step_num > polyploid.general_sim_config.stop_at_step:
        return

    #collect sequences right before WGD, and get them ready to evolve through the next set of trees
    log.write_to_log("\n\n{0}. Collect sequences right before WGD".format(polyploid.analysis_step_num))
    root_seq_files_written_by_gene_tree_by_child_tree = root_seq_maker.run_root_seq_maker(polyploid,
                                             polyploid_subgenomes,
                                             gene_data_by_gt_name, evolver_results_by_gene_tree)
    if polyploid.analysis_step_num > polyploid.general_sim_config.stop_at_step:
        return

    # For the second leg of the autopolyploid sim, we evolve sequences from
    # the time of WGD (say, 300 MYA) to the present, so
    subtrees=["Left","Right"]
    second_leg_random_seeds = [43,99]
    pooled_gene_tree_results_by_tree={}
    pooled_evolver_results_by_tree_by_replicate = {}
    for i in range(0,len(subtrees)):

        subtree=subtrees[i]
        polyploid.subtree_subfolder=subtree
        polyploid_genome_of_interest = polyploid_subgenomes[i]
        log.write_to_log("Processing subgenome " + polyploid_genome_of_interest)

        log.write_to_log("\n\n{0}. Custom GBD model".format(polyploid.analysis_step_num))
        second_leg_gene_tree_results_by_tree_name = custom_GBD_model.run_run_custom_GBD_model_with_root(
            polyploid, species_tree[1+i],postWGD_simulation_leg,
            root_seq_files_written_by_gene_tree_by_child_tree)
        if polyploid.analysis_step_num > polyploid.general_sim_config.stop_at_step:
            return

        log.write_to_log(
            "\n\n{0}. Prune WGD genes that will be dead before the end of the sim".format(polyploid.analysis_step_num))
        second_leg_gene_tree_results_after_pruning_by_gt_name = gene_shedder.shed_genes_only_for_one_branch(
                polyploid, second_leg_gene_tree_results_by_tree_name, polyploid_genome_of_interest, "P1")
        if polyploid.analysis_step_num > polyploid.general_sim_config.stop_at_step:
            return

        log.write_to_log("\n\n{0}. Evolve sequences through gene trees (Evolver)".format(polyploid.analysis_step_num))
        second_leg_evolver_results_by_gene_tree_by_replicate = gene_evolver.run_evolver_with_root_seq(
                polyploid, second_leg_gene_tree_results_after_pruning_by_gt_name, root_seq_files_written_by_gene_tree_by_child_tree,
                second_leg_random_seeds[i])

        pooled_gene_tree_results_by_tree[subtree]=second_leg_gene_tree_results_after_pruning_by_gt_name
        pooled_evolver_results_by_tree_by_replicate[subtree]= second_leg_evolver_results_by_gene_tree_by_replicate


    polyploid.subtree_subfolder = ""

    #note here, there is only one replicate per evolver run
    log.write_to_log("\n\n{0}. Get Ks for trees (Codeml)".format(polyploid.analysis_step_num))
    codeml_results_by_replicate_num = ks_calculator.run_codeml_on_pooled_results(polyploid,
                                                            pooled_gene_tree_results_by_tree,
                                                            pooled_evolver_results_by_tree_by_replicate)

    genomes_of_interest_by_species = {polyploid.species_name: polyploid_subgenomes}
    Ks_results_by_species_by_replicate_num = {polyploid.species_name: codeml_results_by_replicate_num}

    #note, we dont have the outgroup sorted out for the autopolyploid yet.
    log.write_to_log("\n\n{0}. Plot histograms (matplotlib)".format(polyploid.analysis_step_num))
    final_result_files_by_replicate = ks_histogramer.run_Ks_histogramer(polyploid,
                                                                        genomes_of_interest_by_species,
                                                                        Ks_results_by_species_by_replicate_num)


    log.write_to_log("\n\n{0}. Collate results".format(polyploid.analysis_step_num))
    results_organizer.collate_results(polyploid, final_result_files_by_replicate)

    log.write_to_log("\n\n" + polyploid.species_name + " complete")