import log
from pipeline_modules import sagephy_GBD_model, ks_histogramer, ks_calculator, sagephy_tree_relaxer, gene_evolver, \
    species_tree_maker, root_seq_maker, gene_shedder, results_organizer, gene_tree_maker


def run_autosim(polyploid):


    #the autopolyploid has two simulation legs, one before speciation/wgd, and one after,
    preWGD_simulation_leg=polyploid.simulation_legs[0]
    postWGD_simulation_leg=polyploid.simulation_legs[1]
    polyploid_genomes_of_interest = ['P1', 'P2']

    log.write_to_log("\n\n{0}. Make species trees (custom code)".format(polyploid.analysis_step_num))
    species_tree = species_tree_maker.make_species_trees(polyploid)
    if polyploid.analysis_step_num > polyploid.general_sim_config.stop_at_step:
        return

    log.write_to_log("\n\n{0}. Make gene trees given abrupt speciation".format(polyploid.analysis_step_num))
    base_gene_tree_newicks_by_tree_name = gene_tree_maker.make_randomized_gene_trees(polyploid,species_tree[0])
    if polyploid.analysis_step_num > polyploid.general_sim_config.stop_at_step:
        return

    log.write_to_log("\n\n{0}. Add GBD model to gene trees (SaGePhy)".format(polyploid.analysis_step_num))
    #first_leg_time_range=(0,polyploid.FULL_time_MYA-polyploid.WGD_time_MYA)
    #species_trees = {"all":species_tree[0]} #for tha autopolyploid, all gene trees are based ont he same species tree.
    gene_tree_results_by_tree_name = sagephy_GBD_model.run_sagephy(polyploid,preWGD_simulation_leg,
                                                                 base_gene_tree_newicks_by_tree_name )
    if polyploid.analysis_step_num > polyploid.general_sim_config.stop_at_step:
        return

    log.write_to_log("\n\n{0}. Relax gene trees (SaGePhy)".format(polyploid.analysis_step_num))
    relaxed_gene_tree_results = gene_tree_relaxer.relax(polyploid, preWGD_simulation_leg,
                                                        gene_tree_results_by_tree_name)
    if polyploid.analysis_step_num > polyploid.general_sim_config.stop_at_step:
        return

    # For the first leg of the autopolyploid sim, we evolve sequences from
    # what ever the start time was (say, 500 MYA) to the time of WGD (say, 300 MYA)
    first_leg_random_seed=137
    log.write_to_log("\n\n{0}. Evolve sequences through gene trees (Evolver)".format(polyploid.analysis_step_num))
    evolver_results_by_gene_tree = gene_evolver.run_evolver(polyploid, relaxed_gene_tree_results,
                                                            first_leg_random_seed)
    if polyploid.analysis_step_num > polyploid.general_sim_config.stop_at_step:
        return

    #collect sequences right before WGD, and get then ready to evolve through the next set of trees
    log.write_to_log("\n\n{0}. Collect sequences right before WGD".format(polyploid.analysis_step_num))
    root_seq_files_written_by_gene_tree_by_child_tree = root_seq_maker.run_root_seq_maker(polyploid,
                                            relaxed_gene_tree_results, evolver_results_by_gene_tree)
    if polyploid.analysis_step_num > polyploid.general_sim_config.stop_at_step:
        return

    # For the second leg of the autopolyploid sim, we evolve sequences from
    # the time of WGD (say, 300 MYA) to the present, so
    subtrees=["Left","Right"]
    second_leg_random_seeds = [43,99]
    pooled_gene_tree_results_by_tree={}
    pooled_relaxed_gene_tree_results_by_tree={}
    pooled_shed_gene_tree_results_by_tree={}
    pooled_evolver_results_by_tree_by_replicate = {}
    for i in range(0,len(subtrees)):

        subtree=subtrees[i]
        polyploid.subtree_subfolder=subtree
        log.write_to_log("\n\n{0}. Add GBD model to gene trees (SaGePhy)".format(polyploid.analysis_step_num))
        gene_tree_results_by_tree_name = sagephy_GBD_model.run_sagephy_with_root_seq(polyploid, postWGD_simulation_leg,
                                                                     species_tree[1+i],
                                                                     root_seq_files_written_by_gene_tree_by_child_tree)

        log.write_to_log("\n\n{0}. Relax gene trees (SaGePhy)".format(polyploid.analysis_step_num))
        relaxed_gene_tree_results = sagephy_tree_relaxer.relax(polyploid,  postWGD_simulation_leg,
                                                            gene_tree_results_by_tree_name)

        log.write_to_log("\n\n{0}. Prune trees. ".format(polyploid.analysis_step_num) +
              "At every time step post WGD, cull a certain percent of what remains. (custom code)")
        gene_trees_after_gene_shedding = gene_shedder.shed_genes(polyploid,relaxed_gene_tree_results)

        #note, root_seq_files_written_by_gene_tree might not be 1-1
        #with the gt names. If some gt split during the first leg, there might
        #be two seq or more coming off a single gt, which is
        #why there is a list of seq in eac h "GeneTree0.RootSeq.rep0.txt" etc.
        log.write_to_log("\n\n{0}. Evolve sequences through gene trees (Evolver)".format(polyploid.analysis_step_num))
        evolver_results_by_gene_tree_by_replicate = gene_evolver.run_evolver_with_root_seq(
            polyploid, gene_trees_after_gene_shedding, root_seq_files_written_by_gene_tree_by_child_tree,
            second_leg_random_seeds[i])

        pooled_gene_tree_results_by_tree[subtree]=gene_tree_results_by_tree_name
        pooled_relaxed_gene_tree_results_by_tree[subtree]=relaxed_gene_tree_results
        pooled_shed_gene_tree_results_by_tree[subtree]= gene_trees_after_gene_shedding 
        pooled_evolver_results_by_tree_by_replicate[subtree]= evolver_results_by_gene_tree_by_replicate


    polyploid.subtree_subfolder = ""


    #note here, there is only one replicate per evolver run
    log.write_to_log("\n\n{0}. Get Ks for trees (Codeml)".format(polyploid.analysis_step_num))
    codeml_results_by_replicate_num = ks_calculator.run_codeml_on_pooled_results(polyploid,
                                                            pooled_relaxed_gene_tree_results_by_tree,
                                                            pooled_evolver_results_by_tree_by_replicate)

    genomes_of_interest_by_species = {polyploid.species_name: polyploid_genomes_of_interest}
    Ks_results_by_species_by_replicate_num = {polyploid.species_name: codeml_results_by_replicate_num}
    #note, we dont have the outgroup sorted out for the autopolyploid yet.

    log.write_to_log("\n\n{0}. Plot histograms (matplotlib)".format(polyploid.analysis_step_num))
    final_result_files_by_replicate = ks_histogramer.run_Ks_histogramer(polyploid,
                                                                        genomes_of_interest_by_species,
                                                                        Ks_results_by_species_by_replicate_num)

    #Ks_results_by_species_by_replicate_num)

    log.write_to_log("\n\n{0}. Collate results".format(polyploid.analysis_step_num))
    results_organizer.collate_results(polyploid, final_result_files_by_replicate)

    log.write_to_log("\n\n" + polyploid.species_name + " complete")