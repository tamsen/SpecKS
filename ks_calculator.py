import os

def evolver_out_to_sequences(evolver_out_file, expected_seq_names):
    sequences = {}
    with open(evolver_out_file, 'r') as f:

        while (True):

            line = f.readline()

            if not line:
                break

            for s in expected_seq_names:
                if s in line:
                    splat = line.split()
                    seq_name = splat[0]
                    seq_data = "".join(splat[1:len(splat)])

                    if not seq_name in sequences:
                        sequences[seq_name] = []

                    sequences[seq_name].append(seq_data)

    return sequences


def sequences_to_codeml_in(sequences, codeml_in_file):
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

def run_codeml():

    template_codeml_ctl_file="input_templates/template.codeml.ctl"
    dummy_evolver_out_file = "/home/tamsen/Data/Ks_Genome_Simulator/Evolver/mc_3Seq_1000codons_10reps.paml"
    dummy_codeml_in_file = "/home/tamsen/Data/Ks_Genome_Simulator/Codeml/3Seq_1000codons.codonalign.fa"
    expected_seq_names = ['G0_1', 'G1_1']
    my_sequences = evolver_out_to_sequences(dummy_evolver_out_file, expected_seq_names)

    print(my_sequences)
    sequence_files = sequences_to_codeml_in(my_sequences, dummy_codeml_in_file)
    codeml_control_files = write_codeml_control_files(template_codeml_ctl_file, sequence_files)

    cmds = []
    for i in range(0, len(codeml_control_files)):
        ctl_file = codeml_control_files[i]
        base_name = os.path.basename(ctl_file.replace(".codonalign", "").replace(".ctl", ""))
        dir_name = os.path.dirname(ctl_file)
        cmds.append("mkdir " + base_name)
        cmds.append("cd " + base_name)
        cmds.append("codeml " + ctl_file)
        cmds.append("cd .. ")

    for cmd in cmds:
        print(cmd)