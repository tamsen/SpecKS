import os
import unittest
import allosim
import autosim
import main

class TestMain(unittest.TestCase):
    def test_run_sim(self):

        print("Current Working Directory:\t" + os.getcwd())
        current_dir=os.getcwd()
        parent_dir = os.path.dirname(current_dir)

        expected_output_folder="./test_out/test_main"
        if not os.path.exists(expected_output_folder):
            os.makedirs(expected_output_folder)

        test_config_file= os.path.join( parent_dir,"sample_configs","short-test-run.xml")
        mock_args = ["dummy.py",test_config_file]
        conf = main.setup(mock_args)

        # Time since WGD: 5,10, 15,50,100,200 MYA. Total tree length 500 MY. Make allo and autopoly examples.
        list_of_polyploids = main.make_polyploids(conf)

        for polyploid in list_of_polyploids:

            if polyploid.is_allo():
                allosim.run_allosim(polyploid)
            else:
                autosim.run_autosim(polyploid)

        self.assertTrue(os.path.exists(conf.output_folder))  # ad

        #check a final histogram file was created for the allopolyploid and the autopolyploid
        expected_allo_hist= os.path.join(conf.output_folder,"Allopolyploid_1",
                                         "7_ks_histograms","replicate2","Polyploid P1andP2_paml_hist_maxKS5_NG.png")

        self.assertTrue(os.path.exists(expected_allo_hist))

        #TDI, add an extra test for the autopolyploid case, once I get it working again..

if __name__ == '__main__':
    unittest.main()
