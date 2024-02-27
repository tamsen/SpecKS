import os
import shutil
import unittest
from matplotlib import pyplot as plt
from pathlib import Path
import process_wrapper
from pipeline_modules import gene_evolver, gene_tree_data, ks_calculator, ks_histogramer

class GeneEvolverTests3(unittest.TestCase):

    def test_run_evolver_on_simple_newick(self):

        seq_names = ["P1", "P2"]
        gt_root="SpeciesTree.test"
        expected_ks=2.0
        total_num_iterations=10
        num_codons = 1000

        #output locations
        cwd = os.getcwd()
        print("Current Working Directory:\t" + cwd)
        test_out = os.path.join(cwd, "test_out")
        if not os.path.exists(test_out):
            os.makedirs(test_out)
        test_name="EVOLVER3_run_evolver_on_simple_newick"
        evolver_test_out_folder = os.path.join(test_out,  test_name)
        if not os.path.exists(evolver_test_out_folder):
            os.makedirs(evolver_test_out_folder)

        #input locations
        par_dir= Path(__file__).parent
        test_data_folder =  os.path.join(par_dir,"gene_tree_evolver_test_data")
        evolver_control_file =  os.path.join(test_data_folder,"evolver_input_example.dat")
        codeml_control_file =  os.path.join(test_data_folder,"codeml_input_example.ctl")

        #parsing the gene tree
        gt_file = os.path.join(test_data_folder, gt_root + ".tree")
        gt_leaf_map = os.path.join(test_data_folder, gt_root + ".leafmap")
        num_replicates_per_gene_tree = 1
        num_codons = 1000
        random_seed = 3
        gene_tree_result = gene_tree_data.gene_tree_result("gt_name", gt_file, gt_leaf_map)
        evolver_tree_length=gene_tree_result.get_tree_length_as_in_PAML()

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

        for i in range(0, total_num_iterations):
            replicates_subfolder=os.path.join(evolver_test_out_folder,"rep"+str(i))
            if not os.path.exists(replicates_subfolder):
                os.makedirs(replicates_subfolder)
            replicate_fa_file = os.path.join(replicates_subfolder, "evolver_test_rep"+str(i)+".fa")
            write_seq_for_codeml_input(replicate_fa_file, i, sequences_by_seq_name)
            control_file = write_codeml_control_file(codeml_control_file, replicate_fa_file)
            cmd = ["codeml",os.path.basename(control_file)]
            process_wrapper.run_and_wait_on_process(cmd, replicates_subfolder)

    def test_run_evolver(self):

        expected_ks=2.0
        #because in "evolver_input_example.dat" the total tree length is 2.4,
        #ie, a total of ks acculated on one branch is 1.2, and 1.2 on the other.
        #So the two genes are 2.4 apart. 0.4 is non-synon mutuations, and 2.0 is synon.
        #so Ks between them is 2.0


        num_reps=10
        par_dir= Path(__file__).parent
        evolver_control_file =  os.path.join(par_dir,"gene_tree_evolver_test_data",
            "evolver_input_example.dat")
        codeml_control_file =  os.path.join(par_dir,"gene_tree_evolver_test_data",
            "codeml_input_example.ctl")

        #run a simple newick through EVOLVER, with very basic settings.
        cwd = os.getcwd()
        print("cwd:\t" + cwd)
        test_out = os.path.join(cwd, "test_out")
        if not os.path.exists(test_out):
            os.makedirs(test_out)
        evolver_test_out_folder = os.path.join(test_out,  "just_evolver")
        # remove the old output folder
        if os.path.exists(evolver_test_out_folder):
            shutil.rmtree(evolver_test_out_folder)
        os.makedirs(evolver_test_out_folder)

        cmd = ["evolver", "6", evolver_control_file]
        print("cmd path:\t" + evolver_test_out_folder)
        process_wrapper.run_and_wait_on_process(cmd, evolver_test_out_folder)


        # make a .fa input file for codeml
        seq_names=["P1","P2"]
        evolver_out_file1 = os.path.join(evolver_test_out_folder, "mc.paml")
        evolver_out_file2 = os.path.join(evolver_test_out_folder, "mc.txt")
        if os.path.exists(evolver_out_file1):
            sequences_by_seq_name = get_seq_from_evolver_ouput_file(evolver_out_file1, seq_names)
        else:
            sequences_by_seq_name = get_seq_from_evolver_ouput_file(evolver_out_file2, seq_names)
        self.assertEqual(len(sequences_by_seq_name['P1']), num_reps)

        #print(sequences_by_seq_name)
        all_ML_ks_results = []
        all_NG_ks_results = []
        for i in range(0,num_reps):
            codeml_folder=os.path.join(evolver_test_out_folder,"rep"+str(i))
            if not os.path.exists(codeml_folder):
                os.makedirs(codeml_folder)
            fa_file = os.path.join(codeml_folder, "evolver_test_rep"+str(i)+".fa")
            write_seq_for_codeml_input(fa_file,i, sequences_by_seq_name)
            codeml_ctl_file=write_codeml_control_file(codeml_control_file, fa_file)
            #print(fa_file)

            #run CODEML, and see that we get the original settings back.
            cmd = ["codeml", codeml_ctl_file]
            process_wrapper.run_and_wait_on_process(cmd, codeml_folder)

            # read Ks result
            ML_out_file = os.path.join(codeml_folder, "2ML.dS")
            NG_out_file = os.path.join(codeml_folder, "2NG.dS")
            ML_ks_result_data = ks_histogramer.get_Ks_from_file(ML_out_file)
            NG_ks_result_data = ks_histogramer.get_Ks_from_file(NG_out_file)
            ML_ks_results =[ks.ks_between_orthologs for ks in ML_ks_result_data ]
            NG_ks_results =[ks.ks_between_orthologs for ks in NG_ks_result_data ]
            all_ML_ks_results = all_ML_ks_results + ML_ks_results
            all_NG_ks_results = all_NG_ks_results + NG_ks_results

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
        print(seq_data)
        return seq_data


def write_seq_for_codeml_input(fa_file, replicate_idx, sequences_by_name):

    sequence_names =  sequences_by_name.keys()
    with open(fa_file, 'w') as f:
        for seq_name in sequence_names:
            f.writelines(">" + seq_name + "\n")
            f.writelines(sequences_by_name[seq_name][replicate_idx] + "\n")


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
