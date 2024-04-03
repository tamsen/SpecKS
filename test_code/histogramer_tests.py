import os
import shutil
import unittest

from pipeline_modules import ks_histogramer


class HistogramerTests(unittest.TestCase):
    def test_get_Ks_from_GeneTree0_file(self):

        cwd = os.getcwd()
        print("Current Working Directory:\t" + cwd)
        test_out = os.path.join(cwd, "test_out")
        if not os.path.exists(test_out):
            os.makedirs(test_out)

        test_data_folder = "histogramer_test_data"
        ks_file = os.path.join(test_data_folder, "GeneTree0_2ML.dS")
        ks_results = ks_histogramer.get_Ks_from_file(ks_file)

        for ks_result in ks_results:
            print("ks_value is\t" + str(ks_result.round_trip_ks_between_orthologs)
                  + " between leaves " + str(ks_result.ortholog_pair))
            print("row idx\t" + str(ks_result.row_index))
            print("col idx\t" + str(ks_result.col_index))

        self.assertEqual(len(ks_results), 1)
        self.assertEqual(ks_results[0].ortholog_pair, ['G2_2', 'G1_1'])
        self.assertEqual(ks_results[0].round_trip_ks_between_orthologs, 4.7313)

    def test_get_Ks_from_GeneTree4_file(self):

        cwd = os.getcwd()
        print("Current Working Directory:\t" + cwd)
        test_out = os.path.join(cwd, "test_out")
        if not os.path.exists(test_out):
            os.makedirs(test_out)

        test_data_folder = "histogramer_test_data"
        ks_file = os.path.join(test_data_folder, "GeneTree4_2ML.dS")
        ks_results = ks_histogramer.get_Ks_from_file(ks_file)

        for ks_result in ks_results:
            print("ks_value is\t" + str(ks_result.round_trip_ks_between_orthologs)
                  + " between leaves " + str(ks_result.ortholog_pair))
            print("row idx\t" + str(ks_result.row_index))
            print("col idx\t" + str(ks_result.col_index))

        self.assertEqual(len(ks_results), 6)
        self.assertEqual(ks_results[0].ortholog_pair, ['G1_1','G0_1'])
        self.assertEqual(ks_results[0].round_trip_ks_between_orthologs, 3.3035)
        self.assertEqual(ks_results[5].ortholog_pair, ['G5_2', 'G3_2'])
        self.assertEqual(ks_results[5].round_trip_ks_between_orthologs, 2.9798)

if __name__ == '__main__':
    unittest.main()
