from pipeline_modules import gene_tree_maker, ks_histogramer, ks_calculator, gene_tree_relaxer, gene_evolver, \
    species_tree_maker, gene_shedder


def run_allosim(polyploid):

    print("\n\n{0}. Make species trees (custom code)".format(polyploid.analysis_step_num))
    species_trees = species_tree_maker.make_species_trees(polyploid)
    if polyploid.analysis_step_num > polyploid.general_sim_config.stop_at_step:
        return

    print("\n\n{0}. Make gene trees (SaGePhy)".format(polyploid.analysis_step_num))
    gene_tree_results_by_tree_name = gene_tree_maker.run_sagephy(polyploid, species_trees[0])
    if polyploid.analysis_step_num > polyploid.general_sim_config.stop_at_step:
        return

    print("\n\n{0}. Relax gene trees (SaGePhy)".format(polyploid.analysis_step_num))
    relaxed_gene_tree_results = gene_tree_relaxer.relax(polyploid, gene_tree_results_by_tree_name)
    if polyploid.analysis_step_num > polyploid.general_sim_config.stop_at_step:
        return

    print("\n\n{0}. Shed genes post WDG. ".format(polyploid.analysis_step_num) +
          "At every time step post WGD, cull a certain percent of what remains. (custom code)")
    gene_trees_after_gene_shedding = gene_shedder.shed_genes(polyploid, relaxed_gene_tree_results)
    if polyploid.analysis_step_num > polyploid.general_sim_config.stop_at_step:
        return

    print("\n\n{0}. Evolve sequences through gene trees (Evolver)".format(polyploid.analysis_step_num))
    evolver_results_by_gene_tree = gene_evolver.run_evolver(polyploid, relaxed_gene_tree_results)
    if polyploid.analysis_step_num > polyploid.general_sim_config.stop_at_step:
        return

    print("\n\n{0}. Get Ks for trees (Codeml)".format(polyploid.analysis_step_num))
    codeml_results_by_replicate_num = ks_calculator.run_codeml(polyploid,
                                                               relaxed_gene_tree_results, evolver_results_by_gene_tree)
    if polyploid.analysis_step_num > polyploid.general_sim_config.stop_at_step:
        return

    print("\n\n{0}. Plot histograms (matplotlib)".format(polyploid.analysis_step_num))
    ks_histogramer.run_Ks_histogramer(polyploid, codeml_results_by_replicate_num, relaxed_gene_tree_results)
    if polyploid.analysis_step_num > polyploid.general_sim_config.stop_at_step:
        return

    print("\n\n" + polyploid.species_name + " complete.\n\n")