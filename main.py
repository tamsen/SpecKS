import species_tree_maker
import gene_tree_maker
import gene_evolver
import tree_pruner
import ks_calculator
import ks_histogramer
import config



def run_sim():

    configuration = config.SpecKS_config()
    out_dir = configuration.output_folder
    print(out_dir)

    print("1. Make species trees (SaGePhy or by hand)")
    #Time since WGD: 5,10, 15,50,100,200 MYA. Total tree length 500 MY. Make allo and autopoly examples.
    #Example tree (O:500, (P1:200, P2:200):300);
    species_tree = species_tree_maker.get_example_allopolyploid_tree(500,300)
    species_tree_maker.plot_allopolyploid_species_tree(out_dir, 500, 300)

    print("2. Make gene trees (SaGePhy)")
    num_gene_trees=10
    gene_tree_maker.write_gene_tree_commands(species_tree,num_gene_trees,out_dir)

    print("3. Evolve sequences through gene trees (Evolver)")
    gene_evolver.write_evolver_commands(out_dir,"gene_tree_file")

    print("4. Prune trees. at every time step, cull a certain percent of what remains")
    tree_pruner.prune_gene_trees()

    print("5. Get Ks for trees (Codeml)")
    #files needed: seq.codonalign.fa, codeml.ctl
    #Make 1000 gene trees per species tree. Tiley adds that "we simulated variation in the size of the gene families by altering the number of genes at the root of each gene tree by allowing the root age to be exponentially distributed with a mean of 600 Ma"
    ks_calculator.run_codeml()

    print("6. Plot histograms (matplotlib, or Ks rates)")
    ks_histogramer.summarize_Ks(out_dir,"my_species_name",5,False,"NG", 'c', 0.01)

if __name__ == '__main__':
    run_sim()

# See PyCharm help at https://www.jetbrains.com/help/pycharm/
