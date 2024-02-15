import os
import common
import shutil
from pathlib import Path
from Bio import Phylo
from io import StringIO
def write_evolver_control_file(template_dat_file,out_dir, random_seed_odd_integer,
                               num_seq, num_codons, num_replcates, tree_length, newick_tree_string):
    lines_to_write = []
    specs_flags = "NUMSEQ NUMCODONS NUMREPLICATES"
    new_evolver_control_file_name = "mc_{0}Seq_{1}codons_{2}reps.dat".format(
        num_seq, num_codons, tree_length)

    new_dat_file = os.path.join(out_dir, new_evolver_control_file_name)

    #Handle an evolver bug we dont seem to need for later versions.
    #Dont seem to need this for EVOLVER in paml version 4.10.7, June 2023
    #*but* do need it for EVOLVER 4.9, March 2015
    s,vn,vd=get_evolver_version_string(out_dir)
    if (vd[0] < 4.0) or ((vd[0] == 4.0)  and (vd[1] < 10.0) ):
        newick_tree_string = work_around_for_evolver_bug(newick_tree_string)
    #newick_tree_string =work_around_for_no_parenthesis_in_newick(newick_tree_string)

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

            if "RANDOM" in line:
                new_line = line.replace("RANDOM", str(random_seed_odd_integer))

            lines_to_write.append(new_line)

            if not line:
                break

    with open(new_dat_file, 'w') as f:

        for line in lines_to_write:
            f.writelines(line)

    return new_dat_file


def work_around_for_evolver_bug(newick_tree_string):
    fixed_newick_tree_string=newick_tree_string
    if newick_tree_string[-2:] != ");":
        fixed_newick_tree_string = "(" + newick_tree_string.replace(";", ");")
    return fixed_newick_tree_string

def work_around_for_no_parenthesis_in_newick(newick_tree_string):

    open_parenthesis = [i for i in range(0, len(newick_tree_string)) if newick_tree_string[i] == "("]
    if len(open_parenthesis) ==0:
        fixed_newick_tree_string = "(" + newick_tree_string.replace(";", ");")
        newick_tree_string = fixed_newick_tree_string
    return newick_tree_string
def write_evolver_commands(out_dir,random_seed_odd_integer,
                           num_replicates,num_codons,tree_length,gene_tree_result):

    par_dir= Path(__file__).parent.parent
    template_evolver_control_file =  os.path.join(par_dir,"paml_input_templates",
        "template.MCcodon.dat")

    test_tree = Phylo.read(StringIO(gene_tree_result.simple_newick), "newick")

    s,vn,vd=get_evolver_version_string(out_dir)
    if (vd[0] < 4.0) or ((vd[0] == 4.0)  and (vd[1] < 10.0) ):
        X1 = Phylo.to_networkx(gene_tree_result.tree)
        X2 = Phylo.to_networkx(test_tree)
        nodes1 = list(X1.nodes)
        nodes2 = list(X2.nodes)
        num_seq = len(nodes1)
        print("nodes1=\t" + str(len(nodes1)))
        print("nodes2=\t" + str(len(nodes2)))
        print("num seq = num nodes:\t" + str(num_seq))
    else:
        #for recent versions of PAML
        num_seq = gene_tree_result.num_terminal_leaves
        print("num seq = num_extant_leaves:\t" + str(num_seq))

    evolver_control_file = write_evolver_control_file(template_evolver_control_file,
                                                      out_dir, random_seed_odd_integer,
                                                      num_seq, num_codons,
                                                      num_replicates, tree_length,
                                                      gene_tree_result.simple_newick)

    cmd= ["evolver","6", evolver_control_file]
    return cmd
def get_evolver_version_string(gene_tree_subfolder):
    cmd = ["evolver","-v"]
    out_string, error_string = common.run_and_wait_on_process(cmd, gene_tree_subfolder)
    version_string=out_string.split("\n")[0]
    version_number=version_string.split()[4][0:-1]
    version_decimals=[int(v) for v in version_number.split(".")]
    print(version_string)
    print("v"+ version_number)
    print("v:" + str(version_decimals))
    return version_string, version_number, version_decimals

def run_evolver_with_root_seq(polyploid, gene_tree_results_by_gene_tree_name,
                              root_seq_files_written_by_gene_tree, tree_length_for_this_leg, random_seed_odd_integer):

    config = polyploid.general_sim_config
    if len(polyploid.subtree_subfolder) > 0:
        subfolder = os.path.join(polyploid.species_subfolder,
                                 str(polyploid.analysis_step_num) + "_sequence_evolver_" + polyploid.subtree_subfolder)
    else:
        subfolder = os.path.join(polyploid.species_subfolder, str(polyploid.analysis_step_num) + "_sequence_evolver")

    print(subfolder)
    if not os.path.exists(subfolder):
        os.makedirs(subfolder)

    evolver_tree_length = get_evolver_tree_length(config, tree_length_for_this_leg)
    evolver_results_by_gene_tree_by_replicate={}
    for gene_tree_name,gene_tree_result in gene_tree_results_by_gene_tree_name.items():

        evolver_results_by_replicate= {}
        gene_tree_subfolder=os.path.join(subfolder,gene_tree_result.gene_tree_name)
        os.makedirs(gene_tree_subfolder)
        sequence_replicates=root_seq_files_written_by_gene_tree[gene_tree_name]

        for replicate_i in range(0,(len(sequence_replicates))):
            replicate_subfolder = os.path.join(gene_tree_subfolder , "rep" + str(replicate_i))
            replicate_seq_file = sequence_replicates[replicate_i]
            os.makedirs(replicate_subfolder)
            dst=os.path.join(replicate_subfolder, "RootSeq.txt")
            shutil.copyfile(replicate_seq_file, dst)

            cmd = write_evolver_commands(replicate_subfolder,random_seed_odd_integer,
                                         1,
                                         config.num_codons, evolver_tree_length, gene_tree_result)
            out_string,error_string = common.run_and_wait_on_process(cmd, replicate_subfolder)

            evolver_result_file_A=os.path.join(replicate_subfolder,"mc.txt")
            evolver_result_file_B=os.path.join(replicate_subfolder,"mc.paml")

            if os.path.exists(evolver_result_file_A):
                evolver_results_by_replicate[replicate_i] =evolver_result_file_A
                #evolver_results_by_gene_tree[gene_tree_name]=evolver_result_file_A
            elif os.path.exists(evolver_result_file_B):
                evolver_results_by_replicate[replicate_i] =evolver_result_file_B
                #evolver_results_by_gene_tree[gene_tree_name]=evolver_result_file_B
            else:
                error_string="Evolver failed to output a sequence file.  It should be in " + \
                              gene_tree_subfolder + " but its not. Check the evolver stderr."
                print("Error: " + error_string)
                raise ValueError(error_string)

        evolver_results_by_gene_tree_by_replicate[gene_tree_name]=evolver_results_by_replicate

    polyploid.analysis_step_num=polyploid.analysis_step_num+1
    return evolver_results_by_gene_tree_by_replicate



def run_evolver(polyploid, gene_tree_results_by_gene_tree_name, tree_length_for_this_leg, random_seed_odd_integer):


    config = polyploid.general_sim_config
    if len(polyploid.subtree_subfolder) > 0:
        subfolder = os.path.join(polyploid.species_subfolder,
                                 str(polyploid.analysis_step_num) + "_sequence_evolver_" + polyploid.subtree_subfolder)
    else:
        subfolder = os.path.join(polyploid.species_subfolder, str(polyploid.analysis_step_num) + "_sequence_evolver")

    print(subfolder)
    if not os.path.exists(subfolder):
        os.makedirs(subfolder)

    evolver_tree_length = get_evolver_tree_length(config, tree_length_for_this_leg)
    evolver_results_by_gene_tree={}
    for gene_tree_name,gene_tree_result in gene_tree_results_by_gene_tree_name.items():

        gene_tree_subfolder=os.path.join(subfolder,gene_tree_result.gene_tree_name)
        os.makedirs(gene_tree_subfolder)

        print("gene tree file:\t " + gene_tree_result.gene_tree_file_name)
        print("\t\tnewick:\t " + gene_tree_result.simple_newick)
        print("\t\tnum leaves:\t " + str(gene_tree_result.num_terminal_leaves))

        if gene_tree_result.num_terminal_leaves < 2:
            #Then all we have left is the outgroup. No point in running evolver.
            continue

        cmd = write_evolver_commands(gene_tree_subfolder, random_seed_odd_integer,
                                     config.num_replicates_per_gene_tree,
                                     config.num_codons, evolver_tree_length, gene_tree_result)
        common.run_and_wait_on_process(cmd, gene_tree_subfolder)

        #common evolver complaints:
        #if "perhaps too many '('?" in error_string:  <- fix newick
        # See this with evolver 4.9. Where is asks for "# seq" it actually wants the num nodes
        # (including internal nodes). Where as more recent versions want the num extant leaves.


        #Error: error in tree: too many species in tree. <- make sure the requested seq matches the tree topology

        #this mess is b/c diff versions of evoler can output diff file names, A & B
        evolver_result_file_A=os.path.join(gene_tree_subfolder,"mc.txt")
        evolver_result_file_B=os.path.join(gene_tree_subfolder,"mc.paml")
        if os.path.exists(evolver_result_file_A):
            evolver_results_by_gene_tree[gene_tree_name]=evolver_result_file_A
        elif os.path.exists(evolver_result_file_B):
            evolver_results_by_gene_tree[gene_tree_name]=evolver_result_file_B
        else:
            error_string="Evolver failed to output a sequence file.  It should be in " + \
                             gene_tree_subfolder + " but its not. Check the evolver stderr."
            print("Error: " + error_string)
            raise ValueError(error_string)

    polyploid.analysis_step_num=polyploid.analysis_step_num+1
    return evolver_results_by_gene_tree


def get_evolver_tree_length(config, tree_length_for_this_leg):
    # The "evolver" tree length is the expected number of substitutions per site along all
    # branches in the phylogeny, calculated as the sum of the branch lengths
    # so, if my tree length is 500MY, and Ks = 0.01/Myr”, =>
    # Then, the "evolver" tree length is 0.01*500 = 5.
    evolver_tree_length = config.Ks_per_Myr * tree_length_for_this_leg
    return evolver_tree_length

