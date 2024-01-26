import os
import common
def write_evolver_control_file(template_dat_file, new_dat_file,
                               num_seq, num_codons, num_replcates, tree_length, newick_tree_string):
    lines_to_write = []
    specs_flags = "NUMSEQ NUMCODONS NUMREPLICATES"

    newick_tree_string = work_around_for_evolver_bug(newick_tree_string)

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
                print("using newick " + newick_tree_string)

            if "LENGTH" in line:
                new_line = line.replace("LENGTH", str(tree_length))

            lines_to_write.append(new_line)

            if not line:
                break

    with open(new_dat_file, 'w') as f:

        for line in lines_to_write:
            f.writelines(line)

    return new_dat_file


def work_around_for_evolver_bug(newick_tree_string):

    if newick_tree_string[-2:] != ");":
        fixed_newick_tree_string = "(" + newick_tree_string.replace(";", ");")
        newick_tree_string = fixed_newick_tree_string
    return newick_tree_string


def write_evolver_commands(out_dir,num_replicates,num_codons,tree_length,gene_tree_result):

    template_evolver_control_file="./input_templates/template.MCcodon.dat"
    num_seq=gene_tree_result.info_dict["No. of vertices"]
    new_evolver_control_file_name="mc_{0}Seq_{1}codons_{2}reps.dat".format(
        num_seq,num_codons,tree_length)

    new_evolver_control_file=os.path.join(out_dir,new_evolver_control_file_name)


    evolver_control_file = write_evolver_control_file(template_evolver_control_file,
                                                      new_evolver_control_file,
                                                      num_seq, num_codons,
                                                      num_replicates, tree_length,
                                                      gene_tree_result.simple_newick)

    cmd= ["evolver","6", evolver_control_file]
    return cmd

def run_evolver(config, gene_tree_results_by_gene_tree_name):

    out_dir = config.output_folder
    subfolder = os.path.join(out_dir, str(config.sim_step_num) + "_sequence_evolver")
    if not os.path.exists(subfolder):
        os.makedirs(subfolder)

    evolver_results_by_gene_tree={}
    for gene_tree_name,gene_tree_result in gene_tree_results_by_gene_tree_name.items():

        gene_tree_subfolder=os.path.join(subfolder,gene_tree_result.gene_tree_name)
        os.makedirs(gene_tree_subfolder)

        print("gene tree file:\t " + gene_tree_result.gene_tree_file_name)
        print("\t\tnewick:\t " + gene_tree_result.simple_newick)
        print("\t\tnum leaves:\t " + str(gene_tree_result.num_extant_leaves))
        cmd = write_evolver_commands(gene_tree_subfolder,config.num_replicates_per_gene_tree,
                                     config.num_codons,config.tree_length, gene_tree_result)
        common.run_and_wait_on_process(cmd, gene_tree_subfolder)

        evolver_result_file_A=os.path.join(gene_tree_subfolder,"mc.txt")
        evolver_result_file_B=os.path.join(gene_tree_subfolder,"mc.paml")
        if os.path.exists(evolver_result_file_A):
            evolver_results_by_gene_tree[gene_tree_name]=evolver_result_file_A
        elif os.path.exists(evolver_result_file_B):
            evolver_results_by_gene_tree[gene_tree_name]=evolver_result_file_B
        else:
            raise ValueError('Evolver failed to output a sequence file.  It should be in '
                             + gene_tree_subfolder + " but its not. Check the evolver stderr.")

    config.sim_step_num=config.sim_step_num+1
    return evolver_results_by_gene_tree

