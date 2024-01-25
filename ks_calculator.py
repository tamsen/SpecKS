import os

import common


def evolver_out_to_sequences(evolver_out_file, expected_seq_names):
    sequences_by_leaf = {}
    with open(evolver_out_file, 'r') as f:

        while (True):

            line = f.readline()

            if not line:
                break

            for s in expected_seq_names:
                if s in line:
                    splat = line.split()
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
    replacement_flags = ["DUMMY.codonalign.fa", "mlcTree_DUMMY.out"]

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


def write_codeml_control_files(template_ctl_file, sequence_files):
    results = [write_codeml_control_file(template_ctl_file, sf) for sf in sequence_files]
    return results

def run_codeml(config,evolver_results_by_gene_tree):

    out_dir = config.output_folder
    subfolder = os.path.join(out_dir, "5_ks_calcs_by_codeml")
    if not os.path.exists(subfolder):
        os.makedirs(subfolder)

    template_codeml_ctl_file="input_templates/template.codeml.ctl"
    codeml_results_by_replicate_num={}
    replicates = [r for r in range(0,config.num_replicates_per_gene_tree)]

    for gene_tree_name, evolver_output_file in evolver_results_by_gene_tree.items():

        gene_tree_subfolder = os.path.join(subfolder, gene_tree_name)
        if not os.path.exists(gene_tree_subfolder):
            os.makedirs(gene_tree_subfolder)

        expected_seq_names = ['G0_1', 'G1_1']

        # sequences_by_leaf is plural bc we did replicates
        sequences_by_leaf = evolver_out_to_sequences(evolver_output_file, expected_seq_names)
        sequence_files_written=[]

        for r in replicates:
            replicates_subfolder = os.path.join(gene_tree_subfolder,"replicate_"+ str(r))
            if not os.path.exists(replicates_subfolder):
                os.makedirs(replicates_subfolder)

            replicate_fa_file = os.path.join(replicates_subfolder,"3Seq_1000codons.codonalign.rep" + str(r) + ".fa")
            #error, this file is coming up empty:
            # 3Seq_1000codons.codonalign.rep4.fa
            #specifically failing for gene tree 1

            with open(replicate_fa_file, 'w') as f:

                for seq_name in expected_seq_names:
                    if seq_name in sequences_by_leaf:
                        f.writelines(">" + seq_name + "\n")
                        f.writelines(sequences_by_leaf[seq_name][r] + "\n")

            sequence_files_written.append(replicate_fa_file )
            control_file = write_codeml_control_file(template_codeml_ctl_file, replicate_fa_file)
            cmd = ["codeml",control_file]
            print(cmd)
            common.run_and_wait_on_process(cmd, replicates_subfolder)
            result= codeml_result(replicates_subfolder)

            if r not in codeml_results_by_replicate_num:
               codeml_results_by_replicate_num[r]={}

            codeml_results_by_replicate_num[r][gene_tree_name]= result

    return codeml_results_by_replicate_num

class codeml_result():
    NG_file=""
    ML_file=""
    def __init__(self, gene_tree_subfolder):
        self.NG_file = os.path.join(gene_tree_subfolder,"2NG.dN")
        self.ML_file = os.path.join(gene_tree_subfolder,"2ML.dN")