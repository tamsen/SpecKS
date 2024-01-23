import os
import subprocess

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

def write_evolver_commands(out_dir,gene_tree_newick):

    template_evolver_control_file="input_templates/template.MCcodon.dat"
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
                                                      num_replicates, tree_length,
                                                      gene_tree_newick)

    cmd= ["evolver","6", evolver_control_file]
    print(cmd)
    return cmd

def run_evolver(config,gene_trees_by_file_name):

    out_dir = config.output_folder
    subfolder = os.path.join(out_dir, "evolver")
    if not os.path.exists(subfolder):
        os.makedirs(subfolder)

    print("my env:")
    print(str(os.environ))
    #my env: 'VIRTUAL_ENV': '/home/tamsen/Git/SpecKS/SpecKS/venv'
    #the env that has what I need:   / opt / anaconda3 / envs / mynewenv / bin / python
    #https: // stackoverflow.com / questions / 7438681 / how - to - duplicate - virtualenv

    gene_tree_files =  list(gene_trees_by_file_name.keys())
    for gene_tree_file in gene_tree_files[0:1]:

        gene_tree_newick = gene_trees_by_file_name[gene_tree_file]
        cmd = write_evolver_commands(out_dir, gene_tree_newick)
        result = subprocess.run(cmd, capture_output=True)
        print("result:\t" + str(result))
        

#need paml in env
#conda activate mynewenv