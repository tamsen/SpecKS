import unittest
from io import StringIO
from Bio import Phylo

from pipeline_modules import custom_GBD_maker

class TestCustomGBD(unittest.TestCase):
    def test_custom_GBD_model(self):
        base_gene_tree_newick='(O:100, (P1:73.36, P2:73.36): 26.64);'
        #custom_GBD_maker.add_GBD_to_GT(base_gene_tree_newick, 1)

        tree = Phylo.read(StringIO(base_gene_tree_newick), "newick")
        clade = tree.clade
        self.iterate_though_clades(clade)

        self.assertEqual(True, False)  # add assertion here

    #if death rate = 0.001359 gene per MY, thats a death rate of 1 per 700 MY
    #
    def get_nodes_to_add(self, clade):
        length=clade.branch_length
        spaces_between_nodes=50
        first_space=spaces_between_nodes/2.0
        node_life_span= 30

        nodes_to_add=[]
        #if length < first_space:

    def iterate_though_clades(self, clade):

        for c in clade.clades:

            if c.name:
                print(c.name +":\t" +str(c.branch_length))
            else:
                print("internal_node" +":\t" +str(c.branch_length))
            self.iterate_though_clades(c)



if __name__ == '__main__':
    unittest.main()
