import os
import shutil
import unittest

from pipeline_modules import ks_histogramer


class HistogramerTests(unittest.TestCase):
    def test_get_Ks_from_file(self):

        cwd = os.getcwd()
        print("Current Working Directory:\t" + cwd)
        test_out = os.path.join(cwd, "test_out")
        if not os.path.exists(test_out):
            os.makedirs(test_out)

        test_data_folder = "histogramer_test_data"
        ks_file = os.path.join(test_data_folder, "GeneTree0_2ML.dS")
        ks_results = ks_histogramer.get_Ks_from_file(ks_file)

        for ks_result in ks_results:
            print("ks_value is\t" + str(ks_result.ks_between_leaves)
                  + " between leaves " + str(ks_result.leaves))
            print("row idx\t" + str(ks_result.row_index))
            print("col idx\t" + str(ks_result.col_index))

        self.assertEqual(True, True)  # add assertion here


if __name__ == '__main__':
    unittest.main()
