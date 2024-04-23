
import os
import unittest
import kp_reader
from kp_classifier import KP_classifier
from kp_data_to_process import data_to_process

class TestClassifyKP(unittest.TestCase):

    def test_ClassifyKP(self):

        known_auto_examples=["TORX","ZUQW","QAUE","TZWR"]
        known_allo_examples=["CSUV"]
        required=["RYJX","UHLI","WAFT","PIVW"]

        data_directory = "/home/tamsen/Data/"
        Ks_file_1 = "final_ks_values_CSUV.fa"
        kp_directory = os.path.join(data_directory, "1KP_final_ks_files")
        out_data_directory = os.path.join(data_directory, "1KP_classifier_out")

        lookup_file = os.path.join(kp_directory, "1kP-Sample-List.csv")
        sample_lookup = kp_reader.get_sample_code_lookup(lookup_file)

        if not os.path.exists(out_data_directory):
            os.makedirs(out_data_directory)

        max_Ks = 2.0
        num_bins = 50
        kernel_size = 2
        classifier = KP_classifier(max_Ks, num_bins, kernel_size,known_auto_examples,known_allo_examples)
        KS_data_files = [Ks_file_1 ] + [f for f in os.listdir(kp_directory) if ".fa" in f]
        max_to_process=0
        i=0
        curated_samples= ['OBUY']#known_auto_examples+known_allo_examples+required

        for Ks_file in KS_data_files:

            Ks_file_path = os.path.join(kp_directory, Ks_file)
            (species_code, species) = kp_reader.get_species_code(Ks_file, sample_lookup)
            #if (i > max_to_process) and not(species_code in curated_samples):
            #    continue

            print("classifying " + str(i) + " of " + str(max_to_process))
            ks_values = kp_reader.get_Ks_from_file(Ks_file_path)
            print("species name: " + species)
            species_data_to_classify = data_to_process(species, species_code, ks_values, max_Ks)

            classifier.classify(species_data_to_classify,out_data_directory)

            i=i+1

        classifier.print_summary_file(out_data_directory)
        classifier.plot_summary_files(out_data_directory)
        self.assertEqual(True, True)


if __name__ == '__main__':
    unittest.main()
