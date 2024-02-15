import os
import unittest

from Bio import Phylo

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
    # The Ks *should* average about 2 if all the parameters are set right
    def test_run_evolver_on_simple_newick(self):

        #newick_string = "(O:500,(P1:300,P2:300):200);"
        cwd=os.getcwd()
        print("Current Working Directory:\t" + cwd)
        test_out= os.path.join(cwd,"test_out")
        if not os.path.exists(test_out):
            os.makedirs(test_out)

        evolver_test_out_folder=os.path.join(test_out,"evolver_test_out")
        if not os.path.exists(evolver_test_out_folder):
            os.makedirs(evolver_test_out_folder)


        test_data_folder="gene_tree_evolver_test_data"
        gt_file=os.path.join(test_data_folder,"GeneTree33.test.tree")
        gt_leaf_map=os.path.join(test_data_folder,"GeneTree33.test.leafmap")
        random_seed_odd_integer=41
        num_replicates_per_gene_tree=1
        num_codons=1000

        gene_tree_result=gene_tree_data.gene_tree_result("gt_name",gt_file,gt_leaf_map)
        tree=gene_tree_result.tree
        X = Phylo.to_networkx(tree)
        nodes = list(X.nodes)
        sum_of_the_branch_lengths=0
        for node in nodes:
            distance=tree.distance(node)
            sum_of_the_branch_lengths =sum_of_the_branch_lengths+distance
            name="no name"
            if node.name:
                name=node.name
            print(name + " distance from root:\t" + str(distance) )

        evolver_tree_length = sum_of_the_branch_lengths * 0.01
        #evolver_tree_length = 6
        all_ML_ks_results=[]
        all_NG_ks_results=[]
        for i in range(0,100):

            random_seed=i*2+1
            ML_ks_results,NG_ks_results = run_evolver_loop(evolver_test_out_folder, evolver_tree_length, gene_tree_result, num_codons,
                              num_replicates_per_gene_tree, random_seed)
            all_ML_ks_results= all_ML_ks_results+ ML_ks_results
            all_NG_ks_results= all_NG_ks_results + NG_ks_results

        print(all_ML_ks_results)
        print(all_NG_ks_results)

        ML_avg=sum(all_ML_ks_results)/len(all_ML_ks_results)
        NG_avg=sum(all_NG_ks_results)/len(all_NG_ks_results)
        print("avg ML " + str(ML_avg))
        print("avg NG " + str(NG_avg))
        #avg ML 1.851786
        #avg NG 1.845688118811881
        self.assertTrue(True)

def run_evolver_loop(evolver_test_out_folder, evolver_tree_length, gene_tree_result, num_codons,
                         num_replicates_per_gene_tree, random_seed_odd_integer):
        # run evolver
        cmd = gene_evolver.write_evolver_commands(evolver_test_out_folder, random_seed_odd_integer,
                                                  num_replicates_per_gene_tree,
                                                  num_codons, evolver_tree_length, gene_tree_result)
        process_wrapper.run_and_wait_on_process(cmd, evolver_test_out_folder)
        # make a .fa input file for codeml
        fa_file = os.path.join(evolver_test_out_folder, "evolver_test.fa")
        evolver_out_file = os.path.join(evolver_test_out_folder, "mc.paml")
        seq_names = ["G1_0", "G2_0"]
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