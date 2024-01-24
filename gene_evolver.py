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

def write_evolver_commands(out_dir,gene_tree_newick):

    cwd=os.getcwd()
    print("cwd:\t" + cwd)

    template_evolver_control_file="./input_templates/template.MCcodon.dat"
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


    gene_tree_files =  list(gene_trees_by_file_name.keys())
    for gene_tree_file in gene_tree_files[0:1]:

        gene_tree_newick = gene_trees_by_file_name[gene_tree_file]
        print("input newick " + gene_tree_newick)
        cmd = write_evolver_commands(subfolder, gene_tree_newick)
        result_CompletedProcess = subprocess.run(cmd, capture_output=True)
        print("cmd:\t" + str(cmd))
        print("result:\t" + str(result_CompletedProcess))
        print("stderr:\t" +result_CompletedProcess.stderr.decode())
        print("stdout:\t" +result_CompletedProcess.stdout.decode())

#current error mesg:
#result:	CompletedProcess(args=['evolver', '6', '/home/tamsen/Data/SpecKS_output/mc_3Seq_1000codons_5reps.dat'], returncode=255,
# stdout=b'EVOLVER in paml version 4.10.7, June 2023\n\n
# Reading options from data file /home/tamsen/Data/SpecKS_output/mc_3Seq_1000codons_5reps.dat\n\nSimulated data will go into mc.txt.\ntranslated aa sequences will go into mc.aa.txt.\n\n3 seqs, 1000 sites, 10 replicate(s)\nSeq file will be about 121K bytes.\n\nSeq #3 (S3) is missing in the tree\n', stderr=b'')


#print("my env:")
#print(str(os.environ))

#need paml in env
#conda activate SpecksEnv1
#open pycharm from the right venv with codeml & evolver to get this to work