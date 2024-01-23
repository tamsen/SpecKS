import species_tree_maker
import gene_tree_maker
import gene_evolver
import tree_pruner
import ks_calculator
import ks_histogramer
import config



def run_sim():

    conf = config.SpecKS_config()
    out_dir = conf.output_folder
    print(out_dir)

    print("1. Make species trees (SaGePhy or by hand)")
    #Time since WGD: 5,10, 15,50,100,200 MYA. Total tree length 500 MY. Make allo and autopoly examples.
    species_tree = species_tree_maker.make_species_trees(conf)

    print("2. Make gene trees (SaGePhy)")
    num_gene_trees=2#10
    gene_trees_by_file_name = gene_tree_maker.run_sagephy(conf, species_tree, num_gene_trees)

    print("3. Evolve sequences through gene trees (Evolver)")
    gene_evolver.run_evolver(conf, gene_trees_by_file_name)

    print("4. Prune trees. at every time step, cull a certain percent of what remains")
    #tree_pruner.prune_gene_trees()

    print("5. Get Ks for trees (Codeml)")
    #files needed: seq.codonalign.fa, codeml.ctl
    #Make 1000 gene trees per species tree. Tiley adds that "we simulated variation in the size of the gene families by altering the number of genes at the root of each gene tree by allowing the root age to be exponentially distributed with a mean of 600 Ma"
    #ks_calculator.run_codeml()

    print("6. Plot histograms (matplotlib, or Ks rates)")
    ks_histogramer.summarize_Ks(out_dir,"my_species_name",5,False,"NG", 'c', 0.01)

if __name__ == '__main__':
    run_sim()

# See PyCharm help at https://www.jetbrains.com/help/pycharm/
