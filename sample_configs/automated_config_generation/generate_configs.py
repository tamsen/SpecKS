import os
import unittest

import process_wrapper
from config import PolyploidParams


class Generate_Config_Files(unittest.TestCase):

    JOBNAME="specks0gsh"
    CONFIGPATH="config=/usr/scratch2/userdata2/tdunn/cluster_commands/sim10_gshed/mesx0-SPC-020.xml"
    def test_making_configs(self):

        template_xml_file="mesx-template.xml"
        template_sh_file="qsub-template.sh"
        local_out_dir="/home/tamsen/Data/SpecKS_input"
        commands_folder_on_mesx="/usr/scratch2/userdata2/tdunn/cluster_commands"
        sim_subfolder="sim11_test"
        destination_folder = os.path.join(commands_folder_on_mesx, sim_subfolder)
        origin_folder=os.path.join(local_out_dir, sim_subfolder)

        if not os.path.exists(local_out_dir):
            os.makedirs(local_out_dir)
        if not os.path.exists(origin_folder):
            os.makedirs(origin_folder)

        decimals_needed=3
        formatter = "{:0" + str(decimals_needed) + "d}"
        spec_times=[80,60,40,20,10]
        wgd_times = [75, 55, 35, 15, 5]
        poly_params_by_name={}

        for i in range(0,len(spec_times)):

            spec_time=spec_times[i]
            wgd_time=wgd_times[i]
            job_number=str(i+1)

            allo_name="Allo"+job_number+"_S"+formatter.format(spec_time)+"W"+formatter.format(wgd_time)
            auto_name="Auto"+job_number+"_S"+formatter.format(spec_time)+"W"+formatter.format(spec_time)

            allo_params = PolyploidParams(spec_time, wgd_time, allo_name)
            auto_params = PolyploidParams(spec_time, spec_time, auto_name)

            poly_params_by_name[allo_name]=allo_params
            poly_params_by_name[auto_name]=auto_params

        for poly_name,poly_params in poly_params_by_name.items():

            new_xml_file_name = poly_name + ".xml"
            xml_replacements=[("POLYPLOID_SECTION",poly_params.to_xml())]
            new_file_created=write_config_file(
                template_xml_file, origin_folder,new_xml_file_name,xml_replacements)
            self.assertTrue(os.path.exists(new_file_created))

            config_path=os.path.join(destination_folder,new_xml_file_name)
            new_script_file_name = poly_name + ".sh"
            sh_replacements=[("JOBNAME",poly_name),("CONFIGPATH",config_path)]
            new_sh_created=write_config_file(
                template_sh_file,origin_folder,new_script_file_name,sh_replacements)
            self.assertTrue(os.path.exists(new_sh_created))

            cmd2=['scp','-r',origin_folder,'tdunn@mesx.sdsu.edu:' +commands_folder_on_mesx]
            print(" ".join(cmd2))
            out_string, error_string = process_wrapper.run_and_wait_on_process(cmd2, local_out_dir)

            #can I qsub remotely?
            #ssh "$server" "mkdir -p $destiny"
            #https://stackoverflow.com/questions/26278167/submitting-a-job-to-qsub-remotely
        #
def write_config_file(template_file, out_dir, new_file_name, replacements):

    lines_to_write = []
    new_config_file = os.path.join(out_dir, new_file_name)

    with open(template_file, 'r') as f:

        while (True):

            line = f.readline()
            new_line = line

            for r_tuple in replacements:

                out_with_the_old=r_tuple[0]
                in_with_the_new=r_tuple[1]

                if out_with_the_old in line:
                    new_line = line.replace(out_with_the_old,in_with_the_new)
                #if "CONFIGPATH" in line:
                #new_line = line.replace("CONFIGPATH", "boo_path")

            lines_to_write.append(new_line)

            if not line:
                break

    with open(new_config_file, 'w') as f:

        for line in lines_to_write:
            f.writelines(line)

    return new_config_file



if __name__ == '__main__':
    unittest.main()
