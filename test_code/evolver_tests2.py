import os
import shutil
import unittest
from matplotlib import pyplot as plt
from pipeline_modules import ks_histogramer
import process_wrapper

class GeneEvolverTests2(unittest.TestCase):

    def test_run_evolver(self):

        expected_ks=2.0
        num_reps=100

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

        evolver_control_file = os.path.join(cwd, "evolver_input_example_2p0_100_bl200_tree.dat")
        cmd = ["evolver", "6", evolver_control_file]
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
        template_ctl_file = os.path.join(cwd, "codeml_input_example3.ctl")
        all_ML_ks_results = []
        all_NG_ks_results = []
        for i in range(0,num_reps):
            codeml_folder=os.path.join(evolver_test_out_folder,"rep"+str(i))
            if not os.path.exists(codeml_folder):
                os.makedirs(codeml_folder)
            fa_file = os.path.join(codeml_folder, "evolver_test_rep"+str(i)+".fa")
            write_seq_for_codeml_input(fa_file,i, sequences_by_seq_name)
            codeml_ctl_file=write_codeml_control_file(template_ctl_file, fa_file)
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
