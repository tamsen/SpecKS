import os
import shutil


def collate_results(polyploid, ks_results_files_by_replicate_num):

    subfolder = os.path.join(polyploid.species_subfolder, str(polyploid.analysis_step_num) + "_final_results")

    if not os.path.exists(subfolder):
        os.makedirs(subfolder)

    for replicate in ks_results_files_by_replicate_num:
        ks_results = ks_results_files_by_replicate_num[replicate]
        for outfile in ks_results:
            dst = os.path.join(subfolder,os.path.basename(outfile))
            shutil.copyfile(outfile, dst)