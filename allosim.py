import gene_evolver
import gene_tree_maker
import gene_tree_relaxer
import ks_calculator
import ks_histogramer
import species_tree_maker
import tree_pruner


def run_allosim(polyploid):

    print("\n\n{0}. Make species trees (custom code)".format(polyploid.analysis_step_num))
    species_trees = species_tree_maker.make_species_trees(polyploid)

    print("\n\n{0}. Make gene trees (SaGePhy)".format(polyploid.analysis_step_num))
    gene_tree_results_by_tree_name = gene_tree_maker.run_sagephy(polyploid, species_trees[0])

    print("\n\n{0}. Relax gene trees (SaGePhy)".format(polyploid.analysis_step_num))
    relaxed_gene_tree_results = gene_tree_relaxer.relax(polyploid, gene_tree_results_by_tree_name)

    print("\n\n{0}. Evolve sequences through gene trees (Evolver)".format(polyploid.analysis_step_num))
    evolver_results_by_gene_tree = gene_evolver.run_evolver(polyploid, relaxed_gene_tree_results)

    print("\n\n{0}. Prune trees. ".format(polyploid.analysis_step_num) +
          "At every time step post WGD, cull a certain percent of what remains. (custom code)")
    tree_pruner.prune_gene_trees(polyploid)

    print("\n\n{0}. Get Ks for trees (Codeml)".format(polyploid.analysis_step_num))
    codeml_results_by_replicate_num = ks_calculator.run_codeml(polyploid,
                                                               relaxed_gene_tree_results, evolver_results_by_gene_tree)

    print("\n\n{0}. Plot histograms (matplotlib)".format(polyploid.analysis_step_num))
    ks_histogramer.run_Ks_histogramer(polyploid, codeml_results_by_replicate_num, relaxed_gene_tree_results)

    print("\n\n" + polyploid.species_name + " complete.\n\n")