# https://stackoverflow.com/questions/58919356/how-to-save-phylo-tree-to-file
# https://biopython.org/wiki/Phylo

from Bio import Phylo
import matplotlib.pyplot as plt
from io import StringIO

def save_tree_plot(newick_string, out_file_name):

    tree = Phylo.read(StringIO(newick_string), "newick")
    fig = plt.figure(figsize=(10, 20), dpi=100)
    axes = fig.add_subplot(1, 1, 1)
    Phylo.draw(tree, axes=axes, do_show=False)
    # plt.show()
    plt.savefig(out_file_name)
    #simple_newick=out_file_name.replace(".png",".simple.tree")
    #Phylo.write([tree], simple_newick, "newick")

    plt.clf()
    plt.close()