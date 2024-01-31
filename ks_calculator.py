import os
import shutil

import common


def evolver_out_to_sequences(evolver_out_file):
    sequences_by_leaf = {}
    with open(evolver_out_file, 'r') as f:

        while (True):

            line = f.readline()

            if not line:
                break

            splat = line.split()

            if len(splat) < 1:
                continue

            if "_" not in line:
                continue

            seq_name_aka_leaf_name = splat[0]
            seq_data = "".join(splat[1:len(splat)])

            if not seq_name_aka_leaf_name in sequences_by_leaf:
                sequences_by_leaf[seq_name_aka_leaf_name] = []

            sequences_by_leaf[seq_name_aka_leaf_name].append(seq_data)

    return sequences_by_leaf




def sequences_to_codeml_in(sequences, codeml_in_file):

    if len(sequences)== 0:
        return []
    sequence_names = list(sequences.keys())
    num_replicates = len(sequences[sequence_names[0]])
    files_written = []

    for i in range(0, num_replicates):

        replicate_out_file = codeml_in_file.replace(".fa", ".rep" + str(i) + ".fa")

        with open(replicate_out_file, 'w') as f:

            for seq_name in sequence_names:
                f.writelines(">" + seq_name + "\n")
                f.writelines(sequences[seq_name][i] + "\n")

        files_written.append(replicate_out_file)

    return files_written

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


def run_codeml_on_pooled_results(polyploid, pooled_relaxed_gene_tree_results,
                                 pooled_evolver_results_by_subtree_by_gene_tree_by_replicate):

    config = polyploid.general_sim_config
    subfolder = os.path.join(polyploid.species_subfolder, str(polyploid.analysis_step_num) + "_ks_calcs_by_codeml")
    if not os.path.exists(subfolder):
        os.makedirs(subfolder)

    template_codeml_ctl_file= "paml_input_templates/template.codeml.ctl"
    codeml_results_by_replicate_num={}
    replicates = [r for r in range(0,config.num_replicates_per_gene_tree)]
    subtree_names = pooled_evolver_results_by_subtree_by_gene_tree_by_replicate.keys()
    gene_tree_names = pooled_evolver_results_by_subtree_by_gene_tree_by_replicate["Right"].keys()

    for r in replicates:
        codeml_results_by_replicate_num[r]={}

    for gene_tree_name in gene_tree_names:

        gene_tree_subfolder = os.path.join(subfolder, gene_tree_name)
        if not os.path.exists(gene_tree_subfolder):
            os.makedirs(gene_tree_subfolder)

            sequence_files_written=[]
            for r in replicates:
                print("\t replicate " + str(r))
                replicates_subfolder = os.path.join(gene_tree_subfolder, "replicate_" + str(r))
                if not os.path.exists(replicates_subfolder):
                    os.makedirs(replicates_subfolder)

                    replicate_fa_file = os.path.join(replicates_subfolder,
                                                     gene_tree_name + ".codons.rep" + str(r) + ".fa")

                    for subtree in subtree_names:
                        dst = os.path.join(replicates_subfolder, subtree + "_mc.paml")
                        pooled_evolver_results_by_subtree_by_gene_tree=pooled_evolver_results_by_subtree_by_gene_tree_by_replicate[subtree][gene_tree_name]

                        if r in pooled_evolver_results_by_subtree_by_gene_tree:
                            evolver_out =pooled_evolver_results_by_subtree_by_gene_tree[r]
                            shutil.copyfile(evolver_out, dst)

                            #TODO - put all these seq together into an .fa file that codeml can take
                            #read the sequnces we want. Note, there should only be one seq per leaf
                            species_of_interest = ['P1', 'P2']
                            sequences_by_leaf = get_sequences_for_leaves_within_the_polyploid(
                                evolver_out, gene_tree_name, pooled_relaxed_gene_tree_results[subtree],
                                species_of_interest)
                        else:
                            continue

                        with open(replicate_fa_file, 'a') as f:
                            for seq_name, seqs in sequences_by_leaf.items():
                                for seq in seqs:
                                    seq_name_with_subtree=seq_name+"_"+subtree
                                    f.writelines(">" + seq_name_with_subtree + "\n")
                                    f.writelines(seq + "\n")

                    sequence_files_written.append(replicate_fa_file )
                    control_file = write_codeml_control_file(template_codeml_ctl_file, replicate_fa_file)
                    cmd = ["codeml",os.path.basename(control_file)]
                    print("\t cmd: " + " ".join(cmd))
                    print("\t cwd: " + replicates_subfolder)
                    print("\t calculating Ks.. ")
                    common.run_and_wait_on_process(cmd, replicates_subfolder)
                    print("\t Ks determined...")
                    result= codeml_result(replicates_subfolder)
                    codeml_results_by_replicate_num[r][gene_tree_name] = result

    return codeml_results_by_replicate_num
def run_codeml(polyploid,relaxed_gene_tree_results, evolver_results_by_gene_tree):

    config = polyploid.general_sim_config
    subfolder = os.path.join(polyploid.species_subfolder, str(polyploid.analysis_step_num) + "_ks_calcs_by_codeml")
    if not os.path.exists(subfolder):
        os.makedirs(subfolder)

    template_codeml_ctl_file= "paml_input_templates/template.codeml.ctl"
    codeml_results_by_replicate_num={}
    replicates = [r for r in range(0,config.num_replicates_per_gene_tree)]

    for gene_tree_name, evolver_output_file in evolver_results_by_gene_tree.items():

        gene_tree_subfolder = os.path.join(subfolder, gene_tree_name)
        if not os.path.exists(gene_tree_subfolder):
            os.makedirs(gene_tree_subfolder)

        print("calculationg Ks for " + gene_tree_name)

        species_of_interest=['P1','P2']
        sequences_by_leaf = get_sequences_for_leaves_within_the_polyploid(
            evolver_output_file, gene_tree_name,relaxed_gene_tree_results,species_of_interest)
        sequence_files_written=[]

        for r in replicates:
            print("\t replicate " + str(r))
            replicates_subfolder = os.path.join(gene_tree_subfolder,"replicate_"+ str(r))
            if not os.path.exists(replicates_subfolder):
                os.makedirs(replicates_subfolder)

            replicate_fa_file = os.path.join(replicates_subfolder,
                                             gene_tree_name+".codons.rep" + str(r) + ".fa")

            with open(replicate_fa_file, 'w') as f:

                for seq_name, seqs in sequences_by_leaf.items():
                        f.writelines(">" + seq_name + "\n")
                        f.writelines(seqs[r] + "\n")

            if len(sequences_by_leaf.keys())==0:
                print("Warning. No gene tree leaves made it into the polyploid for htis replicate.")
                continue
 
            sequence_files_written.append(replicate_fa_file )
            control_file = write_codeml_control_file(template_codeml_ctl_file, replicate_fa_file)
            cmd = ["codeml",os.path.basename(control_file)]
            print("\t cmd: " + " ".join(cmd))
            print("\t cwd: " + replicates_subfolder)
            print("\t calculating Ks.. ")
            common.run_and_wait_on_process(cmd, replicates_subfolder)
            print("\t Ks determined...")
            result= codeml_result(replicates_subfolder)

            if r not in codeml_results_by_replicate_num:
               codeml_results_by_replicate_num[r]={}

            codeml_results_by_replicate_num[r][gene_tree_name]= result

    polyploid.analysis_step_num=polyploid.analysis_step_num+1
    return codeml_results_by_replicate_num


def get_sequences_for_leaves_within_the_polyploid(
        evolver_output_file, gene_tree_name, relaxed_gene_tree_results,species_of_interest):

    # sequences_by_leaf is plural bc we did replicates
    sequences_by_leaf = evolver_out_to_sequences(evolver_output_file)

    #Get the map of which leaves in the gene tree are in which species
    #gene_tree_name='GeneTree0'
    #but pooled_relaxed_gene_tree_results have "_Left" or "_Right" appended to each key name
    leaf_map = relaxed_gene_tree_results[gene_tree_name].leaves_by_species

    #Get only leaves which are in the polyploid (ie P1 and P2)
    leafs_of_interest = []
    for species, leaves_in_species in leaf_map.items():
        if species in species_of_interest: #["P1", "P2"]:
            leafs_of_interest = leafs_of_interest + leaves_in_species

    #Now, get the sequences associated with each leaf in the polyploid  (ie P1 and P2).
    #note, there are replicates, so there is a list of sequences for each leaf
    sequences_by_leaf_in_polyploid_only = {}
    for leaf, seq in sequences_by_leaf.items():
        if leaf in leafs_of_interest:
            sequences_by_leaf_in_polyploid_only[leaf] = sequences_by_leaf[leaf]

    sequences_by_leaf = sequences_by_leaf_in_polyploid_only
    return sequences_by_leaf


class codeml_result():
    NG_file=""
    ML_file=""
    def __init__(self, gene_tree_subfolder):
        self.NG_file = os.path.join(gene_tree_subfolder,"2NG.dS")
        self.ML_file = os.path.join(gene_tree_subfolder,"2ML.dS")