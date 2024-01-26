import unittest
import gene_evolver

class GeneEvolverTests(unittest.TestCase):
    def test_newick_fixer(self):

        evolver_OK_1="((G0_0:500.0,(G1_1:200.0,(G2_2:150.095,G3_2:150.09)G4_2:49.90)G5_3:300.0)G6_4:0.0);"
        evolver_OK_2="(O:500, (P1:200, P2:200): 300);"
        evolver_bug="(G0_0:500.0,(G1_1:200.0,(G2_2:150.095,G3_2:150.09)G4_2:49.90)G5_3:300.0)G6_4:0.0;"

        # Note, evolver_bug causes evolver to compain "Error: expecting ; in the tree file."
        # Even though, clearly ";" is there. Adding a few parenthesis seems to fix it..

        fixed_1=gene_evolver.work_around_for_evolver_bug(evolver_OK_1)
        self.assertEqual(fixed_1, evolver_OK_1)

        fixed_2=gene_evolver.work_around_for_evolver_bug(evolver_OK_2)
        self.assertEqual(fixed_2, evolver_OK_2)

        fixed_3=gene_evolver.work_around_for_evolver_bug(evolver_bug)
        self.assertEqual(fixed_3, evolver_OK_1)