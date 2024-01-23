from make_gene_trees import write_SaGePhy_GuestTreeGen_commands


def run_sim(name):

    print("1. Make species trees (SaGePhy or by hand)")
    #Time since WGD: 5,10, 15,50,100,200 MYA. Total tree length 500 MY. Make allo and autopoly examples.
    #Example tree (O:500, (P1:200, P2:200):300);

    print("2. Make gene trees (SaGePhy)")
    #GuestTreeGen [options] [species tree (string or file)] [dup rate] [loss rate] [trans rate] [out prefix] need a species tree from (1), [dup rate],[loss rate] ,[trans rate =0] Dup rate and loss rate should be pulled from a distribution (Tiley: b(4, 2460) and b(4,2053)).

    species_tree = '(O:500,(P1:200,P2:200):300);'
    dup_rate = 0.00162
    loss_rate = 0.00194
    out_file_name = "foo"
    write_SaGePhy_GuestTreeGen_commands(species_tree, dup_rate, loss_rate, out_file_name)

    print("3. Evolve sequences through gene trees (Evolver)")
    #Need: random number seed (odd number), # seqs (leaves on the gene tree), # codons (length of gene being simulated, say 1000 codons), # replicates (~1000). And transition/transversion rate ratio (kappa = 2). And ω = dN/dS rate ratio (omega = 0.2). ANd gene tree topology from the last step.

    print("4. Prune trees. at every time step, cull a certain percent of what remains")

    print("5. Get Ks for trees (Codeml)")
    #files needed: seq.codonalign.fa, codeml.ctl
    #Make 1000 gene trees per species tree. Tiley adds that "we simulated variation in the size of the gene families by altering the number of genes at the root of each gene tree by allowing the root age to be exponentially distributed with a mean of 600 Ma"

    print("6. Plot histograms (matplotlib, or Ks rates)")


if __name__ == '__main__':
    run_sim('Hi SpecKS')

# See PyCharm help at https://www.jetbrains.com/help/pycharm/
