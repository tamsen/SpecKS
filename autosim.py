import gene_evolver
import gene_tree_maker
import gene_tree_relaxer
import ks_calculator
import root_seq_maker
import species_tree_maker
import tree_pruner
import ks_histogramer

def run_autosim(polyploid):

    print("Autopolyploid sim, TODO, in progress...")

    print("\n\n{0}. Make species trees (custom code)".format(polyploid.analysis_step_num))
    species_tree = species_tree_maker.make_species_trees(polyploid)

    print("\n\n{0}. Make gene trees up to WGD (SaGePhy)".format(polyploid.analysis_step_num))
    gene_tree_results_by_tree_name = gene_tree_maker.run_sagephy(polyploid, species_tree[0])

    print("\n\n{0}. Relax gene trees (SaGePhy)".format(polyploid.analysis_step_num))
    relaxed_gene_tree_results = gene_tree_relaxer.relax(polyploid, gene_tree_results_by_tree_name)

    #TODO probably have to set tree length here. It defaults to500
    print("\n\n{0}. Evolve sequences through gene trees (Evolver)".format(polyploid.analysis_step_num))
    evolver_results_by_gene_tree = gene_evolver.run_evolver(polyploid, relaxed_gene_tree_results)

    #collect sequences right before WGD, and get then ready to evolve through the next set of trees
    print("\n\n{0}. Collect sequences right before WGD".format(polyploid.analysis_step_num))
    root_seq_files_written_by_gene_tree = root_seq_maker.run_root_seq_maker(polyploid,relaxed_gene_tree_results, evolver_results_by_gene_tree)

    subtrees=["Left","Right"]
    pooled_gene_tree_results_by_tree={}
    pooled_relaxed_gene_tree_results_by_tree={}
    pooled_evolver_results_by_tree_by_replicate = {}
    for i in range(0,len(subtrees)):

        subtree=subtrees[i]
        polyploid.subtree_subfolder=subtree
        print("\n\n{0}. Make gene trees after WGD (SaGePhy)".format(polyploid.analysis_step_num))
        gene_tree_results_by_tree_name = gene_tree_maker.run_sagephy(polyploid,
                                                                     species_tree[1+i])

        print("\n\n{0}. Relax gene trees (SaGePhy)".format(polyploid.analysis_step_num))
        relaxed_gene_tree_results = gene_tree_relaxer.relax(polyploid, gene_tree_results_by_tree_name)

        # TODO (1) probably have to set tree length here. It defaults to500
        # TODO (2) moove all the root_seq files into the evolver cwd
        print("\n\n{0}. Evolve sequences through gene trees (Evolver)".format(polyploid.analysis_step_num))
        evolver_results_by_gene_tree_by_replicate = gene_evolver.run_evolver_with_root_seq(
            polyploid, relaxed_gene_tree_results, root_seq_files_written_by_gene_tree)

        pooled_gene_tree_results_by_tree[subtree]=gene_tree_results_by_tree_name
        pooled_relaxed_gene_tree_results_by_tree[subtree]=relaxed_gene_tree_results
        pooled_evolver_results_by_tree_by_replicate[subtree]= evolver_results_by_gene_tree_by_replicate


    polyploid.subtree_subfolder = ""

    print("\n\n{0}. Prune trees. ".format(polyploid.analysis_step_num) +
          "At every time step post WGD, cull a certain percent of what remains. (custom code)")
    tree_pruner.prune_gene_trees(polyploid)

    #note here, there is only one replicate per evolver run
    print("\n\n{0}. Get Ks for trees (Codeml)".format(polyploid.analysis_step_num))
    codeml_results_by_replicate_num = ks_calculator.run_codeml_on_pooled_results(polyploid,
                                                               pooled_relaxed_gene_tree_results_by_tree,
                                                               pooled_evolver_results_by_tree_by_replicate)

    #todo - this part not all working yet...
    print("\n\n{0}. Plot histograms (matplotlib)".format(polyploid.analysis_step_num))
    print(codeml_results_by_replicate_num)
    ks_histogramer.run_Ks_histogramer(polyploid, codeml_results_by_replicate_num,
                                      pooled_relaxed_gene_tree_results_by_tree)

    print("\n\n" + polyploid.species_name + " complete")