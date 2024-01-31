import os
import unittest
from pipeline_modules import gene_evolver


class GeneEvolverTests(unittest.TestCase):
    def test_newick_fixer(self):

        evolver_OK_1="((G0_0:500.0,(G1_1:200.0,(G2_2:150.095,G3_2:150.09)G4_2:49.90)G5_3:300.0)G6_4:0.0);"
        evolver_OK_2="(O:500, (P1:200, P2:200): 300);"
        evolver_bug="(G0_0:500.0,(G1_1:200.0,(G2_2:150.095,G3_2:150.09)G4_2:49.90)G5_3:300.0)G6_4:0.0;"

        # Note, evolver_bug causes evolver to compain "Error: expecting ; in the tree file."
        # Even though, clearly ";" is there. Adding a few parenthesis seems to fix it..
        # This issue is not in all evolve versions. Ie, its not in  4.10.7, June 2023

        fixed_1= gene_evolver.work_around_for_evolver_bug(evolver_OK_1)
        self.assertEqual(fixed_1, evolver_OK_1)

        fixed_2= gene_evolver.work_around_for_evolver_bug(evolver_OK_2)
        self.assertEqual(fixed_2, evolver_OK_2)

        fixed_3= gene_evolver.work_around_for_evolver_bug(evolver_bug)
        self.assertEqual(fixed_3, evolver_OK_1)

    def test_check_evolver_version(self):

        test_out="./test_out"
        if not os.path.exists(test_out):
            os.makedirs(test_out)

        version_string, version_number, version_decimals = gene_evolver.get_evolver_version_string(test_out)
        self.assertTrue("version" in version_string)