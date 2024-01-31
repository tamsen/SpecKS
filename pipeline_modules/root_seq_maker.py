import os

from pipeline_modules.ks_calculator import get_sequences_for_leaves_within_the_polyploid


def sequences_to_root_seq_in(sequences, codeml_in_file):

    if len(sequences)== 0:
        return []
    sequence_names = list(sequences.keys())
    num_replicates = len(sequences[sequence_names[0]])
    files_written = []

    for i in range(0, num_replicates):

        replicate_out_file = codeml_in_file.replace(".txt", ".rep" + str(i) + ".txt")

        with open(replicate_out_file, 'w') as f:

            for seq_name in sequence_names:
                f.writelines(sequences[seq_name][i] + "\n")

        files_written.append(replicate_out_file)

    return files_written


def run_root_seq_maker(polyploid, relaxed_gene_tree_results, evolver_results_by_gene_tree):

    config = polyploid.general_sim_config
    subfolder=os.path.join(polyploid.species_subfolder, str(polyploid.analysis_step_num) + "_wgd_root_sequences")
    if not os.path.exists(subfolder):
        os.makedirs(subfolder)

    sequence_files_written_by_gene_tree = {}
    for gene_tree_name, evolver_output_file in evolver_results_by_gene_tree.items():

        gene_tree_subfolder = os.path.join(subfolder, gene_tree_name)
        if not os.path.exists(gene_tree_subfolder):
            os.makedirs(gene_tree_subfolder)

        print("sorting out post-WGD root sequences for " + gene_tree_name)

        species_of_interest=["P"]
        sequences_by_leaf = get_sequences_for_leaves_within_the_polyploid(evolver_output_file, gene_tree_name,
                                                                       relaxed_gene_tree_results,species_of_interest)

        replicate_seq_root_file = os.path.join(gene_tree_subfolder, gene_tree_name + ".RootSeq.txt")
        files_written = sequences_to_root_seq_in(sequences_by_leaf, replicate_seq_root_file)
        sequence_files_written_by_gene_tree[gene_tree_name] = files_written

    polyploid.analysis_step_num = polyploid.analysis_step_num + 1
    print(str(sequence_files_written_by_gene_tree))
    return sequence_files_written_by_gene_tree
