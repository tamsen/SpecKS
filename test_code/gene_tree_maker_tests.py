import unittest

from pipeline_modules import gene_tree_maker


class GeneTreeMakerTests(unittest.TestCase):
    def test_something(self):

        num_gene_trees_needed = 5
        gt_index_formatter = gene_tree_maker.get_gt_index_format(num_gene_trees_needed)
        gt_names=[]
        for i in range(0, num_gene_trees_needed):
            out_file_name = "GeneTree" + gt_index_formatter.format(i)
            gt_names.append(out_file_name)
        self.assertEqual(gt_names[1], "GeneTree1")  # add assertion here
        self.assertEqual(gt_names[4], "GeneTree4")  # add

        num_gene_trees_needed = 1000
        gt_index_formatter = gene_tree_maker.get_gt_index_format(num_gene_trees_needed)
        gt_names=[]

        for i in range(0, num_gene_trees_needed):
            out_file_name = "GeneTree" + gt_index_formatter.format(i)
            gt_names.append(out_file_name)
        self.assertEqual(gt_names[1], "GeneTree0001")  # add assertion here
        self.assertEqual(gt_names[4], "GeneTree0004")  # a
        self.assertEqual(gt_names[999], "GeneTree0999")  # a



if __name__ == '__main__':
    unittest.main()
