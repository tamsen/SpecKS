#import os
#
#import log
#from pipeline_modules.ks_calculator import get_sequences_for_leaves_within_the_polyploid
#
#
# def sequences_to_root_seq_in(sequences, gene_tree_subfolder, gene_tree_name):
#
#     if len(sequences)== 0:
#         return []
#     sequence_names = list(sequences.keys())
#     num_replicates = len(sequences[sequence_names[0]])
#
#     child_tree_seq = {}
#
#     for seq_name in sequence_names:
#
#         files_written = []
#         subtree_path = os.path.join(gene_tree_subfolder, seq_name)
#         if not os.path.exists(subtree_path):
#             os.makedirs(subtree_path)
#
#         for i in range(0, num_replicates):
#
#             out_file_name = gene_tree_name + "_child" + seq_name +".rep" + str(i) + ".txt"
#             replicate_out_file = os.path.join(subtree_path, out_file_name)
#
#             with open(replicate_out_file, 'w') as f:
#                 f.writelines(sequences[seq_name][i] + "\n")
#
#             files_written.append(replicate_out_file)
#         child_tree_seq[seq_name] = files_written
#
#     return child_tree_seq
#
# "Some people wanted to specify the sequence at the root rather than letting the program generate a \
# random sequence. This can be achieved by putting a sequence in the file RootSeq.txt. The sequence \
# cannot have ambiguities or gaps or stop codons. In almost all simulations, it is simply wrong to fix \
# the root sequence, so you should resist the temptation of making the mistake. If you want the \
# simulation to reflect your particular gene, you may estimate parameters under a model from that \
# gene and then simulate data sets using the parameter estimates. - PAML manual "
# def run_root_seq_maker(polyploid, polyploid_genomes_of_interest, relaxed_gene_tree_results, evolver_results_by_gene_tree):
#
#     config = polyploid.general_sim_config
#     subfolder=os.path.join(polyploid.species_subfolder, str(polyploid.analysis_step_num) + "_wgd_root_sequences")
#     if not os.path.exists(subfolder):
#         os.makedirs(subfolder)
#
#     sequence_files_written_by_gene_tree_by_child_tree = {}
#     for gene_tree_name, evolver_output_file in evolver_results_by_gene_tree.items():
#
#         gene_tree_subfolder = os.path.join(subfolder, gene_tree_name)
#         if not os.path.exists(gene_tree_subfolder):
#             os.makedirs(gene_tree_subfolder)
#
#         log.write_to_log("sorting out post-WGD root sequences for " + gene_tree_name)
#
#         species_of_interest=["P"]
#         sequences_by_leaf = get_sequences_for_leaves_within_the_polyploid(evolver_output_file, gene_tree_name,
#                                                                        relaxed_gene_tree_results,species_of_interest)
#
#         files_written_by_child_tree = sequences_to_root_seq_in(sequences_by_leaf, gene_tree_subfolder, gene_tree_name )
#         sequence_files_written_by_gene_tree_by_child_tree[gene_tree_name] = files_written_by_child_tree
#
#     polyploid.analysis_step_num = polyploid.analysis_step_num + 1
#     #print(str(sequence_files_written_by_gene_tree_by_child_tree))
#     return sequence_files_written_by_gene_tree_by_child_tree
