import os
import shutil
import unittest

from Bio import Phylo
from matplotlib import pyplot as plt

import process_wrapper
from pipeline_modules import gene_evolver, gene_tree_data, ks_calculator, ks_histogramer


class GeneEvolverTests(unittest.TestCase):
    def test_newick_fixer(self):

        evolver_OK_1="((G0_0:500.0,(G1_1:200.0,(G2_2:150.095,G3_2:150.09)G4_2:49.90)G5_3:300.0)G6_4:0.0);"
        evolver_OK_2="(O:500, (P1:200, P2:200): 300);"
        evolver_bug="(G0_0:500.0,(G1_1:200.0,(G2_2:150.095,G3_2:150.09)G4_2:49.90)G5_3:300.0)G6_4:0.0;"

        # Note, evolver_bug causes evolver to compain "Error: expecting ; in the tree file."
        # Even though, clearly ";" is there. Adding a few parenthesis seems to fix it..
        # This issue is not in all evolver versions. Ie, its not in  4.10.7, June 2023

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

    # Idea of test: Run evolver along a known length (say 200), then compute the Ks.
    # The Ks *should* average about 2 (or 4 full path between the leaves)
    # if all the parameters are set right
    def test_run_evolver_on_simple_newick(self):

        seq_names = ["P1", "P2"]
        gt_root="SpeciesTree.test"
        expected_ks=2.0*2.0
        total_num_iterations=100

        self.test_evolver_and_ks_results_are_as_expected(expected_ks, gt_root, seq_names, total_num_iterations,
                                                                    "evolver_simple_newick")

    def test_run_evolver_on_more_complex_newick(self):

        seq_names = ["G1_0", "G2_0"]
        gt_root="GeneTree33.test"
        expected_ks=233.79*2.0*0.01
        total_num_iterations=100

        self.test_evolver_and_ks_results_are_as_expected(expected_ks, gt_root, seq_names,
                                                         total_num_iterations, "evolver_complex_newick")


    def test_evolver_and_ks_results_are_as_expected(self, expected_ks,
                                                    gt_root, seq_names, total_num_iterations, test_name):
        # newick_string = "(O:500,(P1:300,P2:300):200);"
        cwd = os.getcwd()
        print("Current Working Directory:\t" + cwd)
        test_out = os.path.join(cwd, "test_out")
        if not os.path.exists(test_out):
            os.makedirs(test_out)
        evolver_test_out_folder = os.path.join(test_out,  test_name)
        # remove the old output folder
        if os.path.exists(evolver_test_out_folder):
            shutil.rmtree(evolver_test_out_folder)
        os.makedirs(evolver_test_out_folder)
        test_data_folder = "gene_tree_evolver_test_data"
        gt_file = os.path.join(test_data_folder, gt_root + ".tree")
        gt_leaf_map = os.path.join(test_data_folder, gt_root + ".leafmap")
        num_replicates_per_gene_tree = 1
        num_codons = 1000
        gene_tree_result = gene_tree_data.gene_tree_result("gt_name", gt_file, gt_leaf_map)
        tree = gene_tree_result.tree
        X = Phylo.to_networkx(tree)
        nodes = list(X.nodes)
        sum_of_the_branch_lengths = 0
        for node in nodes:
            distance = tree.distance(node)
            sum_of_the_branch_lengths = sum_of_the_branch_lengths + distance
            name = "no name"
            if node.name:
                name = node.name
            print(name + " distance from root:\t" + str(distance))
        total_branch_length = tree.total_branch_length()
        print("total_branch_length:\t" + str(total_branch_length))
        evolver_tree_length = total_branch_length * 0.01 * 1.2
        all_ML_ks_results = []
        all_NG_ks_results = []
        for i in range(0, total_num_iterations):
            random_seed = i * 2 + 1
            ML_ks_results, NG_ks_results = run_evolver_loop(evolver_test_out_folder,
                                                            evolver_tree_length, gene_tree_result, num_codons,
                                                            num_replicates_per_gene_tree, random_seed, seq_names)
            all_ML_ks_results = all_ML_ks_results + ML_ks_results
            all_NG_ks_results = all_NG_ks_results + NG_ks_results
        print(all_ML_ks_results)
        print(all_NG_ks_results)
        ML_histogram_out_file_name = os.path.join(evolver_test_out_folder,
                                                  "ML_ks_results" + "_Ks_histogram.png")
        NG_histogram_out_file_name = os.path.join(evolver_test_out_folder,
                                                  "NG_ks_results" + "_Ks_histogram.png")
        histogram_ks_data(all_ML_ks_results, expected_ks, ML_histogram_out_file_name)
        histogram_ks_data(all_NG_ks_results, expected_ks, NG_histogram_out_file_name)
        ML_avg = sum(all_ML_ks_results) / len(all_ML_ks_results)
        NG_avg = sum(all_NG_ks_results) / len(all_NG_ks_results)
        print("avg ML " + str(ML_avg))
        print("avg NG " + str(NG_avg))
        self.assertTrue(os.path.exists(ML_histogram_out_file_name))
        self.assertTrue(os.path.exists(NG_histogram_out_file_name))


def run_evolver_loop(evolver_test_out_folder, evolver_tree_length, gene_tree_result, num_codons,
                         num_replicates_per_gene_tree, random_seed_odd_integer,seq_names):
        # run evolver
        cmd = gene_evolver.write_evolver_commands(evolver_test_out_folder, random_seed_odd_integer,
                                                  num_replicates_per_gene_tree,
                                                  num_codons, evolver_tree_length, gene_tree_result)
        process_wrapper.run_and_wait_on_process(cmd, evolver_test_out_folder)
        # make a .fa input file for codeml
        fa_file = os.path.join(evolver_test_out_folder, "evolver_test.fa")
        evolver_out_file = os.path.join(evolver_test_out_folder, "mc.paml")
        sequences_by_seq_name = get_seq_from_evolver_ouput_file(evolver_out_file, seq_names)
        write_seq_for_codeml_input(fa_file, sequences_by_seq_name)
        # run codeml
        template_codeml_ctl_file = ks_calculator.get_codeml_ctl_template()
        control_file = ks_calculator.write_codeml_control_file(template_codeml_ctl_file, fa_file)
        cmd = ["codeml", os.path.basename(control_file)]
        process_wrapper.run_and_wait_on_process(cmd, evolver_test_out_folder)
        # read Ks result
        ML_out_file = os.path.join(evolver_test_out_folder, "2ML.dS")
        NG_out_file = os.path.join(evolver_test_out_folder, "2NG.dS")
        ML_ks_results = ks_histogramer.get_Ks_from_file(ML_out_file)
        NG_ks_results = ks_histogramer.get_Ks_from_file(NG_out_file)
        print(ML_ks_results)
        print(NG_ks_results)
        return ML_ks_results,NG_ks_results

def get_seq_from_evolver_ouput_file(evolver_out_file, seq_names):
        seq_data = {}
        with open(evolver_out_file, "r") as f:
            lines = f.readlines()
            for l in lines:
                for seq_name in seq_names:
                    if seq_name in l:
                        data = l.split()
                        codons = data[1:len(data)]
                        dna_seq = "".join(codons)
                        seq_data[seq_name] = dna_seq
        print(seq_data)
        return seq_data


def write_seq_for_codeml_input(fa_file, sequences_by_name):

    sequence_names =  sequences_by_name.keys()
    with open(fa_file, 'w') as f:
        for seq_name in sequence_names:
            f.writelines(">" + seq_name + "\n")
            f.writelines(sequences_by_name[seq_name] + "\n")


def histogram_ks_data(ks_data,expected_ks,dist_histogram_out_file_name):

    ks_avg = sum(ks_data) / len(ks_data)

    fig = plt.figure(figsize=(10, 10), dpi=100)
    n_bins=50
    n, bins, patches = plt.hist(ks_data, bins=n_bins, facecolor='b', alpha=0.25, label='histogram data')
    plt.axvline(x=ks_avg, color='r', linestyle='--', label="empirical_avg")
    plt.axvline(x=expected_ks, color='b', linestyle='--', label="expected_ks")
    plt.legend()
    plt.xlabel("Ks")
    plt.ylabel("Num data points")
    plt.title("Ks data histogram")
    plt.savefig(dist_histogram_out_file_name)
    plt.clf()
    plt.close()

    return dist_histogram_out_file_name