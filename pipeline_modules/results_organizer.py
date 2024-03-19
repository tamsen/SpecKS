import os
import shutil

import log


def collate_results(polyploid, ks_results_files_by_species_by_replicate_num):

    subfolder = os.path.join(polyploid.species_subfolder, str(polyploid.analysis_step_num) + "_final_results")
    log.write_to_log("organizing final results...")

    if not os.path.exists(subfolder):
        os.makedirs(subfolder)

    for species, ks_results_files_by_replicate_num in ks_results_files_by_species_by_replicate_num.items():

        for replicate in ks_results_files_by_replicate_num:
            ks_results = ks_results_files_by_replicate_num[replicate]
            for outfile in ks_results:
                file_name= species +"_"+ os.path.basename(outfile)
                dst = os.path.join(subfolder,file_name)
                shutil.copyfile(outfile, dst)