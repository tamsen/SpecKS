from pipeline_modules import gene_tree_maker, ks_histogramer, ks_calculator, gene_tree_relaxer, gene_evolver, \
    species_tree_maker, root_seq_maker, gene_shedder


def run_autosim(polyploid):


    #the autopolyploid has two simulation legs, one before speciation/wgd, and one after,
    preWGD_simulation_leg=polyploid.simulation_legs[0]
    postWGD_simulation_leg=polyploid.simulation_legs[1]

    print("\n\n{0}. Make species trees (custom code)".format(polyploid.analysis_step_num))
    species_tree = species_tree_maker.make_species_trees(polyploid)
    if polyploid.analysis_step_num > polyploid.general_sim_config.stop_at_step:
        return

    print("\n\n{0}. Make gene trees up to WGD (SaGePhy)".format(polyploid.analysis_step_num))
    first_leg_time_range=(0,polyploid.FULL_time_MYA-polyploid.WGD_time_MYA)
    gene_tree_results_by_tree_name = gene_tree_maker.run_sagephy(polyploid,preWGD_simulation_leg,
                                                                 species_tree[0])
    if polyploid.analysis_step_num > polyploid.general_sim_config.stop_at_step:
        return

    print("\n\n{0}. Relax gene trees (SaGePhy)".format(polyploid.analysis_step_num))
    relaxed_gene_tree_results = gene_tree_relaxer.relax(polyploid, preWGD_simulation_leg,
                                                        gene_tree_results_by_tree_name)
    if polyploid.analysis_step_num > polyploid.general_sim_config.stop_at_step:
        return

    # For the first leg of the autopolyploid sim, we evolve sequences from
    # what ever the start time was (say, 500 MYA) to the time of WGD (say, 300 MYA)
    length_of_time_of_first_leg_of_sim=polyploid.FULL_time_MYA-polyploid.WGD_time_MYA
    print("\n\n{0}. Evolve sequences through gene trees (Evolver)".format(polyploid.analysis_step_num))
    evolver_results_by_gene_tree = gene_evolver.run_evolver(polyploid, relaxed_gene_tree_results,
                                                            length_of_time_of_first_leg_of_sim)
    if polyploid.analysis_step_num > polyploid.general_sim_config.stop_at_step:
        return

    #collect sequences right before WGD, and get then ready to evolve through the next set of trees
    print("\n\n{0}. Collect sequences right before WGD".format(polyploid.analysis_step_num))
    root_seq_files_written_by_gene_tree = root_seq_maker.run_root_seq_maker(polyploid,
                                            relaxed_gene_tree_results, evolver_results_by_gene_tree)
    if polyploid.analysis_step_num > polyploid.general_sim_config.stop_at_step:
        return

    # For the second leg of the autopolyploid sim, we evolve sequences from
    # the time of WGD (say, 300 MYA) to the present, so
    second_leg_of_sim_time=polyploid.WGD_time_MYA
    second_leg_time_range=(polyploid.FULL_time_MYA-polyploid.WGD_time_MYA,polyploid.FULL_time_MYA)
    subtrees=["Left","Right"]
    pooled_gene_tree_results_by_tree={}
    pooled_relaxed_gene_tree_results_by_tree={}
    pooled_shed_gene_tree_results_by_tree={}
    pooled_evolver_results_by_tree_by_replicate = {}
    for i in range(0,len(subtrees)):

        subtree=subtrees[i]
        polyploid.subtree_subfolder=subtree
        print("\n\n{0}. Make gene trees after WGD (SaGePhy)".format(polyploid.analysis_step_num))
        gene_tree_results_by_tree_name = gene_tree_maker.run_sagephy(polyploid, postWGD_simulation_leg,
                                                                     species_tree[1+i])

        print("\n\n{0}. Relax gene trees (SaGePhy)".format(polyploid.analysis_step_num))
        relaxed_gene_tree_results = gene_tree_relaxer.relax(polyploid,  postWGD_simulation_leg,
                                                            gene_tree_results_by_tree_name)

        print("\n\n{0}. Prune trees. ".format(polyploid.analysis_step_num) +
              "At every time step post WGD, cull a certain percent of what remains. (custom code)")
        gene_trees_after_gene_shedding = gene_shedder.shed_genes(polyploid,relaxed_gene_tree_results)

        print("\n\n{0}. Evolve sequences through gene trees (Evolver)".format(polyploid.analysis_step_num))
        evolver_results_by_gene_tree_by_replicate = gene_evolver.run_evolver_with_root_seq(
            polyploid, gene_trees_after_gene_shedding, root_seq_files_written_by_gene_tree,
            second_leg_of_sim_time)

        pooled_gene_tree_results_by_tree[subtree]=gene_tree_results_by_tree_name
        pooled_relaxed_gene_tree_results_by_tree[subtree]=relaxed_gene_tree_results
        pooled_shed_gene_tree_results_by_tree[subtree]= gene_trees_after_gene_shedding 
        pooled_evolver_results_by_tree_by_replicate[subtree]= evolver_results_by_gene_tree_by_replicate


    polyploid.subtree_subfolder = ""


    #note here, there is only one replicate per evolver run
    print("\n\n{0}. Get Ks for trees (Codeml)".format(polyploid.analysis_step_num))
    codeml_results_by_replicate_num = ks_calculator.run_codeml_on_pooled_results(polyploid,
                                                            pooled_relaxed_gene_tree_results_by_tree,
                                                            pooled_evolver_results_by_tree_by_replicate)

    print("\n\n{0}. Plot histograms (matplotlib)".format(polyploid.analysis_step_num))
    print(codeml_results_by_replicate_num)
    ks_histogramer.run_Ks_histogramer(polyploid, codeml_results_by_replicate_num,
                                      pooled_relaxed_gene_tree_results_by_tree)

    print("\n\n" + polyploid.species_name + " complete")