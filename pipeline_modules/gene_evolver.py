import os
import common
import shutil
from pathlib import Path
from Bio import Phylo
from io import StringIO
def write_evolver_control_file(template_dat_file,out_dir,
                               num_seq, num_codons, num_replcates, tree_length, newick_tree_string):
    lines_to_write = []
    specs_flags = "NUMSEQ NUMCODONS NUMREPLICATES"
    new_evolver_control_file_name = "mc_{0}Seq_{1}codons_{2}reps.dat".format(
        num_seq, num_codons, tree_length)

    new_dat_file = os.path.join(out_dir, new_evolver_control_file_name)

    #Handle an evolver bug we dont seem to need for later versions.
    #Dont seem to need this for EVOLVER in paml version 4.10.7, June 2023
    s,vn,vd=get_evolver_version_string(out_dir)
    if (vd[0] < 4.0) or ((vd[0] == 4.0)  and (vd[1] < 10.0) ):
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

    par_dir= Path(__file__).parent.parent
    template_evolver_control_file =  os.path.join(par_dir,"paml_input_templates",
        "template.MCcodon.dat")


    #s,vn,vd=get_evolver_version_string(out_dir)
    #if (vd[0] < 4.0) or ((vd[0] == 4.0)  and (vd[1] < 10.0) ):
    #    num_seq=gene_tree_result.info_dict["No. of vertices"]
    #else:
    #    num_seq = gene_tree_result.info_dict["No. of extant leaves"]

    #todo - carry Tree around with result, instead of just  newick string
    newick = gene_tree_result.simple_newick
    tree= Phylo.read(StringIO(newick), "newick")
    terminal_leaf_names = [t.name for t in  tree.get_terminals()]
    num_seq = len(terminal_leaf_names)

    evolver_control_file = write_evolver_control_file(template_evolver_control_file,
                                                      out_dir,
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
                              root_seq_files_written_by_gene_tree, tree_length_for_this_leg):

    config = polyploid.general_sim_config
    if len(polyploid.subtree_subfolder) > 0:
        subfolder = os.path.join(polyploid.species_subfolder,
                                 str(polyploid.analysis_step_num) + "_sequence_evolver_" + polyploid.subtree_subfolder)
    else:
        subfolder = os.path.join(polyploid.species_subfolder, str(polyploid.analysis_step_num) + "_sequence_evolver")

    print(subfolder)
    if not os.path.exists(subfolder):
        os.makedirs(subfolder)

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

            cmd = write_evolver_commands(replicate_subfolder, 1,
                                         config.num_codons, tree_length_for_this_leg, gene_tree_result)
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



def run_evolver(polyploid, gene_tree_results_by_gene_tree_name, tree_length):

    config = polyploid.general_sim_config
    if len(polyploid.subtree_subfolder) > 0:
        subfolder = os.path.join(polyploid.species_subfolder,
                                 str(polyploid.analysis_step_num) + "_sequence_evolver_" + polyploid.subtree_subfolder)
    else:
        subfolder = os.path.join(polyploid.species_subfolder, str(polyploid.analysis_step_num) + "_sequence_evolver")

    print(subfolder)
    if not os.path.exists(subfolder):
        os.makedirs(subfolder)

    evolver_results_by_gene_tree={}
    for gene_tree_name,gene_tree_result in gene_tree_results_by_gene_tree_name.items():

        gene_tree_subfolder=os.path.join(subfolder,gene_tree_result.gene_tree_name)
        os.makedirs(gene_tree_subfolder)

        print("gene tree file:\t " + gene_tree_result.gene_tree_file_name)
        print("\t\tnewick:\t " + gene_tree_result.simple_newick)
        print("\t\tnum leaves:\t " + str(gene_tree_result.num_terminal_leaves))
        cmd = write_evolver_commands(gene_tree_subfolder,config.num_replicates_per_gene_tree,
                                     config.num_codons,tree_length, gene_tree_result)
        common.run_and_wait_on_process(cmd, gene_tree_subfolder)

        #common evolver complaints:
        #if "perhaps too many '('?" in error_string:  <- fix newick
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

