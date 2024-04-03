import os

import log
import process_wrapper
import shutil
from pathlib import Path
from Bio import Phylo
from io import StringIO

from visualization import tree_visuals_by_phylo


def write_evolver_control_file(template_dat_file, out_dir, evolver_version_is_old, random_seed_odd_integer,
                               num_seq, num_codons, num_replicates, tree_length, newick_tree_string):
    lines_to_write = []
    specs_flags = "NUMSEQ NUMCODONS NUMREPLICATES"
    new_evolver_control_file_name = "mc_{0}Seq_{1}codons_{2}replicates_{3}treelength.dat".format(
        num_seq, num_codons, num_replicates, tree_length)

    new_dat_file = os.path.join(out_dir, new_evolver_control_file_name)

    # Handle an evolver bug we dont seem to need for later versions.
    # Dont seem to need this for EVOLVER in paml version 4.10.7, June 2023
    # *but* do need it for EVOLVER 4.9, March 2015
    if evolver_version_is_old:
        newick_tree_string = work_around_for_evolver_bug(newick_tree_string)

    with open(template_dat_file, 'r') as f:

        while (True):

            line = f.readline()
            new_line = line

            if specs_flags in line:
                new_line = line.replace("NUMSEQ", str(num_seq))
                new_line = new_line.replace("NUMCODONS", str(num_codons))
                new_line = new_line.replace("NUMREPLICATES", str(num_replicates))

            if "TREE" in line:
                new_line = line.replace("TREE", newick_tree_string)
                log.write_to_log("using newick " + newick_tree_string)

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
    fixed_newick_tree_string = newick_tree_string
    if newick_tree_string[-2:] != ");":
        fixed_newick_tree_string = "(" + newick_tree_string.replace(";", ");")
    return fixed_newick_tree_string

def write_evolver_commands(out_dir, template_evolver_control_file, evolver_version_is_old, random_seed_odd_integer,
                           num_replicates, num_codons, tree_length, gene_tree_result):

    # test_tree = Phylo.read(StringIO(gene_tree_result.simple_newick), "newick")

    if evolver_version_is_old:
        num_seq = len(gene_tree_result.get_named_nodes())
    else:
        num_seq = gene_tree_result.num_terminal_leaves

    evolver_control_file = write_evolver_control_file(template_evolver_control_file,
                                                      out_dir,evolver_version_is_old,
                                                      random_seed_odd_integer,
                                                      num_seq, num_codons,
                                                      num_replicates, tree_length,
                                                      gene_tree_result.simple_newick)

    cmd = ["evolver", "6", evolver_control_file]
    return cmd


def get_evolver_version_string(process_folder):
    cmd = ["evolver", "-v"]
    out_string, error_string = process_wrapper.run_and_wait_on_process(cmd, process_folder)
    version_string = out_string.split("\n")[0]
    version_number = version_string.split()[4][0:-1]
    version_decimals = [int(v) for v in version_number.split(".")]
    return evolver_version_string(version_string, version_number, version_decimals)


def run_evolver_with_root_seq(polyploid, gene_tree_results_by_gene_tree_name,
                              root_seq_files_written_by_gene_tree_by_child_tree,
                              random_seed_odd_integer):
    config = polyploid.general_sim_config
    if len(polyploid.subtree_subfolder) > 0:
        subfolder = os.path.join(polyploid.species_subfolder,
                                 str(polyploid.analysis_step_num) + "_sequence_evolver_" + polyploid.subtree_subfolder)
    else:
        subfolder = os.path.join(polyploid.species_subfolder, str(polyploid.analysis_step_num) + "_sequence_evolver")

    if not os.path.exists(subfolder):
        os.makedirs(subfolder)

    template_evolver_control_file = get_evolver_template_file()
    evolver_results_by_gene_tree_by_replicate = {}
    evolver_version_is_old = get_evolver_version_string(subfolder).is_old_version()

    for full_gene_tree_name, gene_tree_result in gene_tree_results_by_gene_tree_name.items():

        splat = full_gene_tree_name.split("_child")
        first_leg_GT_name = splat[0]
        second_leg_GT_name = splat[1]
        evolver_results_by_replicate = {}
        gene_tree_subfolder = os.path.join(subfolder, gene_tree_result.gene_tree_name)
        os.makedirs(gene_tree_subfolder)
        root_seq_files_written_for_replicates = root_seq_files_written_by_gene_tree_by_child_tree[first_leg_GT_name][
            second_leg_GT_name]
        evolver_tree_length = get_evolver_tree_length(config, gene_tree_result)


        for replicate_i in range(0, (len(root_seq_files_written_for_replicates))):
            replicate_subfolder = os.path.join(gene_tree_subfolder, "rep" + str(replicate_i))
            replicate_seq_file = root_seq_files_written_for_replicates[replicate_i]
            os.makedirs(replicate_subfolder)
            dst = os.path.join(replicate_subfolder, "RootSeq.txt")
            shutil.copyfile(replicate_seq_file, dst)

            cmd = write_evolver_commands(replicate_subfolder, template_evolver_control_file,
                                         evolver_version_is_old,
                                         random_seed_odd_integer,
                                         1,
                                         config.num_codons, evolver_tree_length, gene_tree_result)
            out_string, error_string = process_wrapper.run_and_wait_on_process(cmd, replicate_subfolder)

            evolver_result_file_A = os.path.join(replicate_subfolder, "mc.txt")
            evolver_result_file_B = os.path.join(replicate_subfolder, "mc.paml")

            if os.path.exists(evolver_result_file_A):
                evolver_results_by_replicate[replicate_i] = evolver_result_file_A
            elif os.path.exists(evolver_result_file_B):
                evolver_results_by_replicate[replicate_i] = evolver_result_file_B
            else:
                error_string = "Evolver failed to output a sequence file.  It should be in " + \
                               gene_tree_subfolder + " but its not. Check the evolver stderr."
                print("Error: " + error_string)
                log.write_to_log("Error: " + error_string)
                raise ValueError(error_string)

            # GeneTree05_childG6_1.rep1.txt
            # combined_gene_tree_name=first_leg_GT_name +"_child"+child_tree
            evolver_results_by_gene_tree_by_replicate[full_gene_tree_name] = evolver_results_by_replicate

    polyploid.analysis_step_num = polyploid.analysis_step_num + 1
    return evolver_results_by_gene_tree_by_replicate


def run_evolver(polyploid, gene_tree_results_by_gene_tree_name, random_seed_odd_integer):

    config = polyploid.general_sim_config
    include_visualizations = config.include_visualizations

    if len(polyploid.subtree_subfolder) > 0:
        subfolder = os.path.join(polyploid.species_subfolder,
                                 str(polyploid.analysis_step_num) + "_sequence_evolver_" + polyploid.subtree_subfolder)
    else:
        subfolder = os.path.join(polyploid.species_subfolder, str(polyploid.analysis_step_num) + "_sequence_evolver")

    if not os.path.exists(subfolder):
        os.makedirs(subfolder)

    evolver_results_by_gene_tree = {}
    random_seed = random_seed_odd_integer
    template_evolver_control_file = get_evolver_template_file()
    evolver_version_is_old = get_evolver_version_string(subfolder).is_old_version()
    for gene_tree_name, gene_tree_result in gene_tree_results_by_gene_tree_name.items():

        random_seed = random_seed + 2
        gene_tree_subfolder = os.path.join(subfolder, gene_tree_result.gene_tree_name)
        os.makedirs(gene_tree_subfolder)
        evolver_tree_length = get_evolver_tree_length(config, gene_tree_result)
        #log.write_to_log("gene tree file:\t " + gene_tree_result.gene_tree_file_name)
        #print("\t\tnewick:\t " + gene_tree_result.simple_newick)
        #print("\t\tnum leaves:\t " + str(gene_tree_result.num_terminal_leaves))
        #print("\t\tevolver tree length:\t " + str(evolver_tree_length))

        if include_visualizations:
            plot_file_name_1 = os.path.join(gene_tree_subfolder, "gt_used_by_evolver_phylo.png")
            tree_visuals_by_phylo.save_tree_plot_from_newick(gene_tree_result.simple_newick, plot_file_name_1)

        if gene_tree_result.num_terminal_leaves < 2:
            # Then all we have left is the outgroup. No point in running evolver.
            continue

        cmd = write_evolver_commands(gene_tree_subfolder, template_evolver_control_file,
                                     evolver_version_is_old,
                                     random_seed,
                                     config.num_replicates_per_gene_tree,
                                     config.num_codons, evolver_tree_length, gene_tree_result)
        process_wrapper.run_and_wait_on_process(cmd, gene_tree_subfolder)

        # common evolver complaints:
        # if "perhaps too many '('?" in error_string:  <- fix newick
        # See this with evolver 4.9. Where is asks for "# seq" it actually wants the num nodes
        # (including internal nodes). Where as more recent versions want the num extant leaves.

        # Error: error in tree: too many species in tree. <- make sure the requested seq matches the tree topology

        # this mess is b/c diff versions of evoler can output diff file names, A & B
        evolver_result_file_A = os.path.join(gene_tree_subfolder, "mc.txt")
        evolver_result_file_B = os.path.join(gene_tree_subfolder, "mc.paml")
        if os.path.exists(evolver_result_file_A):
            evolver_results_by_gene_tree[gene_tree_name] = evolver_result_file_A
        elif os.path.exists(evolver_result_file_B):
            evolver_results_by_gene_tree[gene_tree_name] = evolver_result_file_B
        else:
            error_string = "Evolver failed to output a sequence file.  It should be in " + \
                           gene_tree_subfolder + " but its not. Check the evolver stderr."
            print("Error: " + error_string)
            raise ValueError(error_string)

    polyploid.analysis_step_num = polyploid.analysis_step_num + 1
    return evolver_results_by_gene_tree


def get_evolver_template_file():
    par_dir = Path(__file__).parent.parent
    template_evolver_control_file = os.path.join(par_dir, "paml_input_templates",
                                                 "evolver_input_example.dat")
    return template_evolver_control_file


def get_evolver_tree_length(config, gene_tree_result):
    # Example calculation:
    # Newick: ((O:500, (P1:100, P2:100)G2_0:400));
    # because in "evolver_input_example.dat" the total tree length is
    # 100+ 100 + 400 + 500 = 1100
    # ie, if we want P1 and P2 to be 1 unit of Ks apart, for 1 MYA
    # Then we expect tree distance of 100+100=200 to have Ks of 1,ie K= Ks+Kn of 1.2
    # So the ratio of tree distance (200) to K (1.2) to be 1.2/200 = 0.006
    # So the K for the whole tree should be (0.006)*1100

    # n general, the ratio of tree distance to K is
    # config.per_site_evolutionary_distance * 1.2 (because K= Ks+Kn , and our Kn/Ks ratio is 0.2) * 0.5 (because we t only one way, not round trip)
    ratio_of_tree_distance_to_K = config.per_site_evolutionary_distance * 0.6

    # (config.Ks_per_Myr  + 0.2 Kn_per_Myr )* 0.5
    # ratio_of_tree_distance_to_K = config.Ks_per_Myr * 0.65

    total_tree_length = gene_tree_result.get_tree_length_as_in_PAML()
    evolver_tree_length = total_tree_length * ratio_of_tree_distance_to_K
    return evolver_tree_length


class evolver_version_string():
    version_string = ""
    version_number = 0
    version_decimals = 0.0

    def __init__(self, vs, vn, vd):
        self.version_string = vs
        self.version_number = vn
        self.version_decimals = vd

    def is_old_version(self):

        log.write_to_log("Evolver version:")
        log.write_to_log(self.version_string)

        if (self.version_decimals[0] < 4.0) or ((self.version_decimals[0] == 4.0)
                                                and (self.version_decimals[1] < 10.0)):
            return True
        else:
            return False
