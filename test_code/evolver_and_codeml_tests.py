import os
import shutil
import unittest
from matplotlib import pyplot as plt
from pathlib import Path
import process_wrapper
from pipeline_modules import gene_evolver, gene_tree_data, ks_calculator, ks_histogramer

class GeneEvolverAndCodemlTests(unittest.TestCase):

    def test_run_evolver_on_simple_newicks_for_a_range_of_Ks(self):

        #self.run_evolver_on_simple_newick("SpeciesTreeKs4.test",4.0,3)
        #self.run_evolver_on_simple_newick("SpeciesTreeKs4.test", 4.0, 1995)
        #self.run_evolver_on_simple_newick("SpeciesTreeKs4.test", 4.0, 43)

        #self.run_evolver_on_simple_newick("SpeciesTreeKs2.test",2.0,3)
        #self.run_evolver_on_simple_newick("SpeciesTreeKs2.test", 2.0, 1995)
        #self.run_evolver_on_simple_newick("SpeciesTreeKs2.test", 2.0, 43)

        self.run_evolver_on_simple_newick("SpeciesTreeKs2_extranode.test",2.0,3)
        #self.run_evolver_on_simple_newick("SpeciesTreeKs2_extranode.test", 2.0, 1995)
        #self.run_evolver_on_simple_newick("SpeciesTreeKs2_extranode.test", 2.0, 43)

        #self.run_evolver_on_simple_newick("SpeciesTreeKs1.test",1.0,3)
        #self.run_evolver_on_simple_newick("SpeciesTreeKs1.test", 1.0, 1995)
        #self.run_evolver_on_simple_newick("SpeciesTreeKs1.test", 1.0, 43)
    def run_evolver_on_simple_newick(self, genetree_name, expected_ks, random_seed):

        seq_names = ["P1", "P2"]
        #gt_root="SpeciesTreeKs2.test"
        #random_seed = 3
        num_replicates_per_gene_tree=100
        num_codons = 1000


        #situation A
        # Newick: ((O:500, (P1:100, P2:100)G2_0:400));
        #because in "evolver_input_example.dat" the total tree length is
        # 1+ 1 + 4 + 5 = 11
        #ie, if we want P1 and P2 to be 2 unit of Ks apart, for 1 MYA
        #Then we expect tree length of 11*1.2 = 13.2
        #so Ks between them is 2.0 (with Ks + Kn being 2.4)

        #situation B
        # Newick: ((O:500, (P1:200, P2:200)G2_0:300));
        # 2 + 2 + 3 + 5 = 12
        #expected_ks=4.0
        #Then we expect tree length of 12*1.2 = 14.4


        #output locations
        cwd = os.getcwd()
        print("Current Working Directory:\t" + cwd)
        test_out = os.path.join(cwd, "test_out")
        if not os.path.exists(test_out):
            os.makedirs(test_out)
        test_name="EVOLVER3_run_evolver_on_" +genetree_name + "randomseed_" + str(random_seed)
        evolver_test_out_folder = os.path.join(test_out,  test_name)
        if not os.path.exists(evolver_test_out_folder):
            os.makedirs(evolver_test_out_folder)

        #input locations
        par_dir= Path(__file__).parent
        test_data_folder =  os.path.join(par_dir,"gene_tree_evolver_test_data")
        evolver_control_file =  os.path.join(test_data_folder,"evolver_input_example.dat")
        codeml_control_file =  os.path.join(test_data_folder,"codeml_input_example.ctl")

        #parsing the gene tree
        gt_file = os.path.join(test_data_folder, genetree_name + ".tree")
        gt_leaf_map = os.path.join(test_data_folder, genetree_name + ".leafmap")
        gene_tree_result = gene_tree_data.gene_tree_result("gt_name", gt_file, gt_leaf_map)
        evolver_tree_length= gene_tree_result.get_tree_length_as_in_PAML()* 0.01 * 1.2

        cmd = gene_evolver.write_evolver_commands(evolver_test_out_folder,
                                                  evolver_control_file,
                                                  random_seed,
                                                  num_replicates_per_gene_tree ,
                                                  num_codons, evolver_tree_length, gene_tree_result)
        process_wrapper.run_and_wait_on_process(cmd, evolver_test_out_folder)
        evolver_out_file1 = os.path.join(evolver_test_out_folder, "mc.paml")
        evolver_out_file2 = os.path.join(evolver_test_out_folder, "mc.txt")
        if os.path.exists(evolver_out_file1):
                sequences_by_seq_name = get_seq_from_evolver_ouput_file(evolver_out_file1, seq_names)
        else:
                sequences_by_seq_name = get_seq_from_evolver_ouput_file(evolver_out_file2, seq_names)

        all_ML_ks_results = []
        all_NG_ks_results = []
        for i in range(0, num_replicates_per_gene_tree):
            replicates_subfolder=os.path.join(evolver_test_out_folder,"rep"+str(i))
            if not os.path.exists(replicates_subfolder):
                os.makedirs(replicates_subfolder)
            replicate_fa_file = os.path.join(replicates_subfolder, "evolver_test_rep"+str(i)+".fa")
            write_seq_for_codeml_input(replicate_fa_file, i, sequences_by_seq_name)
            control_file = write_codeml_control_file(codeml_control_file, replicate_fa_file)
            cmd = ["codeml",os.path.basename(control_file)]
            process_wrapper.run_and_wait_on_process(cmd, replicates_subfolder)

            ML_out_file = os.path.join(replicates_subfolder, "2ML.dS")
            NG_out_file = os.path.join(replicates_subfolder, "2NG.dS")
            ML_ks_result_data = ks_histogramer.get_Ks_from_file(ML_out_file)
            NG_ks_result_data = ks_histogramer.get_Ks_from_file(NG_out_file)
            ML_ks_results = [ks.ks_between_orthologs for ks in ML_ks_result_data]
            NG_ks_results = [ks.ks_between_orthologs for ks in NG_ks_result_data]

            all_ML_ks_results = all_ML_ks_results + ML_ks_results
            all_NG_ks_results = all_NG_ks_results + NG_ks_results

        print(all_ML_ks_results)
        print(all_NG_ks_results)

        xmax=expected_ks*1.5
        ML_histogram_out_file_name = os.path.join(evolver_test_out_folder,
                                                  "ML_ks_results_rand" + str(random_seed)+
                                                  "_xmax" +str(xmax) +
                                                  "_Ks_histogram.png")
        NG_histogram_out_file_name = os.path.join(evolver_test_out_folder,
                                                  "NG_ks_results_rand" + str(random_seed)+
                                                  "_xmax" + str(xmax) +
                                                  "_Ks_histogram.png")

        histogram_ks_data(all_ML_ks_results, expected_ks, ML_histogram_out_file_name,xmax)
        histogram_ks_data(all_NG_ks_results, expected_ks, NG_histogram_out_file_name,xmax)

        ML_histogram_out_file_name = os.path.join(evolver_test_out_folder,
                                                  "ML_ks_results_rand" + str(random_seed) +
                                                  "_Ks_histogram.png")
        NG_histogram_out_file_name = os.path.join(evolver_test_out_folder,
                                                  "NG_ks_results_rand" + str(random_seed) +
                                                  "_Ks_histogram.png")

        histogram_ks_data(all_ML_ks_results, expected_ks, ML_histogram_out_file_name)
        histogram_ks_data(all_NG_ks_results, expected_ks, NG_histogram_out_file_name)

        ML_avg = sum(all_ML_ks_results) / len(all_ML_ks_results)
        NG_avg = sum(all_NG_ks_results) / len(all_NG_ks_results)
        print("avg ML " + str(ML_avg))
        print("avg NG " + str(NG_avg))
        self.assertTrue(os.path.exists(ML_histogram_out_file_name))
        self.assertTrue(os.path.exists(NG_histogram_out_file_name))


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

                        if not seq_name in seq_data:
                            seq_data[seq_name] = []
                        seq_data[seq_name].append(dna_seq)
        #print(seq_data)
        return seq_data


def write_seq_for_codeml_input(fa_file, replicate_idx, sequences_by_name):

    sequence_names =  sequences_by_name.keys()
    with open(fa_file, 'w') as f:
        for seq_name in sequence_names:
            f.writelines(">" + seq_name + "\n")
            f.writelines(sequences_by_name[seq_name][replicate_idx] + "\n")


def histogram_ks_data(ks_data,expected_ks,dist_histogram_out_file_name, xmax=False):

    ks_avg = sum(ks_data) / len(ks_data)

    fig = plt.figure(figsize=(10, 10), dpi=100)
    n_bins=50
    n, bins, patches = plt.hist(ks_data, bins=n_bins, facecolor='b', alpha=0.25, label='histogram data')
    plt.axvline(x=ks_avg, color='r', linestyle='--', label="empirical_avg")
    plt.axvline(x=expected_ks, color='b', linestyle='--', label="expected_ks")
    plt.legend()

    if xmax:
        plt.xlim(0, xmax)
    plt.xlabel("Ks")
    plt.ylabel("Num data points")
    plt.title("Ks data histogram")
    plt.savefig(dist_histogram_out_file_name)
    plt.clf()
    plt.close()

    return dist_histogram_out_file_name

def write_codeml_control_file(template_ctl_file, sequence_file):
    lines_to_write = []
    base_name = os.path.basename(sequence_file).replace(".codonalign", "").replace(".fa", "")
    new_ctl_file = sequence_file.replace(".fa", ".ctl")

    with open(template_ctl_file, 'r') as f:

        while (True):

            line = f.readline()
            new_line = line

            if "DUMMY.codonalign.fa" in line:
                new_line = line.replace("DUMMY.codonalign.fa", sequence_file)

            if "mlcTree_DUMMY.out" in line:
                new_line = line.replace("DUMMY", base_name)

            lines_to_write.append(new_line)

            if not line:
                break

    with open(new_ctl_file, 'w') as f:

        for line in lines_to_write:
            f.writelines(line)

    return new_ctl_file
