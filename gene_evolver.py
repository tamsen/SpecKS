import os
import common
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

def write_evolver_commands(out_dir,num_replicates,num_codons,tree_length,gene_tree_result):

    template_evolver_control_file="./input_templates/template.MCcodon.dat"
    num_seq=gene_tree_result.num_extant_leaves
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

def run_evolver(config, gene_tree_results_by_file_name):

    out_dir = config.output_folder
    subfolder = os.path.join(out_dir, "3_sequence_evolver")
    if not os.path.exists(subfolder):
        os.makedirs(subfolder)

    evolver_results_by_gene_tree={}
    gene_tree_files =  list(gene_tree_results_by_file_name.keys())
    for gene_tree_file in gene_tree_files[0:2]:

        gene_tree_result = gene_tree_results_by_file_name[gene_tree_file]
        gene_tree_subfolder=os.path.join(subfolder,gene_tree_result.gene_tree_name)
        os.makedirs(gene_tree_subfolder)

        print("gene tree file:\t " + gene_tree_file)
        print("\t\tnewick:\t " + gene_tree_result.simple_newick)
        print("\t\tnum leaves:\t " + str(gene_tree_result.num_extant_leaves))
        cmd = write_evolver_commands(gene_tree_subfolder,config.num_replicates_per_gene_tree,
                                     config.num_codons,config.tree_length, gene_tree_result)
        common.run_and_wait_on_process(cmd, gene_tree_subfolder)

        evolver_result_file=os.path.join(gene_tree_subfolder,"mc.txt")
        evolver_results_by_gene_tree[gene_tree_result.gene_tree_name]=evolver_result_file

    return evolver_results_by_gene_tree

