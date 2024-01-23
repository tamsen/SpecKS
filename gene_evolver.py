import os
import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import glob
from scipy.stats import beta

def write_evolver_control_file(template_dat_file, new_dat_file,
                               num_seq, num_codons, num_replcates, tree_length, newick_tree_string):
    lines_to_write = []
    specs_flags = "NUMSEQ NUMCODONS NUMREPLICATES"

    with open(template_dat_file, 'r') as f:

        while (True):

            line = f.readline()
            new_line = line

            if specs_flags in line:
                new_line = line.replace("NUMSEQ", str(num_seq))
                new_line = new_line.replace("NUMCODONS", str(num_codons))
                new_line = new_line.replace("NUMREPLICATES", str(num_replcates))

            if "TREE" in line:
                new_line = line.replace("TREE", newick_tree_string)

            if "LENGTH" in line:
                new_line = line.replace("LENGTH", str(tree_length))

            lines_to_write.append(new_line)

            if not line:
                break

    with open(new_dat_file, 'w') as f:

        for line in lines_to_write:
            f.writelines(line)

    return new_dat_file

def get_gene_tree_from_file(gene_tree_file):
    return '(O:500,(P1:200,P2:200):300);'

def write_evolver_commands(out_dir,gene_tree_file):

    template_evolver_control_file="input_templates/template.MCcodon.dat"
    gene_tree = get_gene_tree_from_file(gene_tree_file)
    num_seq=3
    num_codons=1000
    num_replicates=10
    tree_length=5
    new_evolver_control_file_name="mc_{0}Seq_{1}codons_{2}reps.dat".format(
        num_seq,num_codons,tree_length)

    new_evolver_control_file=os.path.join(out_dir,new_evolver_control_file_name)


    evolver_control_file = write_evolver_control_file(template_evolver_control_file,
                                                      new_evolver_control_file,
                                                      num_seq, num_codons,
                                                      num_replicates, tree_length, gene_tree)

    cmd= "evolver 6 " + evolver_control_file
    print(cmd)