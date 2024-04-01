import os
import unittest

import process_wrapper
from config import PolyploidParams


class Generate_Config_Files(unittest.TestCase):



    def test_making_configs(self):

        sim_subfolder="sim16_N0p1"
        me_at_remote_URL='tdunn@mesx.sdsu.edu'
        template_xml_file="mesx-template.xml"
        template_sh_file="qsub-template.sh"
        local_out_dir="/home/tamsen/Data/SpecKS_input"
        commands_folder_on_mesx="/usr/scratch2/userdata2/tdunn/cluster_commands"
        output_folder_on_mesx="/usr/scratch2/userdata2/tdunn/SpecKS_Output"
        script_destination_folder_on_mesx = os.path.join(commands_folder_on_mesx, sim_subfolder)
        specks_output_path_on_mesx = os.path.join(output_folder_on_mesx, sim_subfolder)
        specks_output_stub_on_mesx = os.path.join(specks_output_path_on_mesx, "specks")
        origin_folder=os.path.join(local_out_dir, sim_subfolder)

        if not os.path.exists(local_out_dir):
            os.makedirs(local_out_dir)
        if not os.path.exists(origin_folder):
            os.makedirs(origin_folder)

        decimals_needed=3
        formatter = "{:0" + str(decimals_needed) + "d}"
        spec_times=[80,60,40,20,10,5]
        wgd_times = [75, 55, 35, 15, 5,1]

        div_log="lognorm,0.5,5.27"
        div_expon_0p1="expon,0,0.1"

        poly_params_by_name={}
        cluster_cmds=[]
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
            xml_replacements=[("POLYPLOID_SECTION",poly_params.to_xml()),
                              ("OUTPUT_ROOT", specks_output_stub_on_mesx),
                              ("DIV_DIST", div_expon_0p1)
                              ]
            new_file_created=write_config_file(
                template_xml_file, origin_folder,new_xml_file_name,xml_replacements)
            self.assertTrue(os.path.exists(new_file_created))
         #   / usr / scratch2 / userdata2 / tdunn / SpecKS_Output / sim9_gbd_in_auto / specks
            config_path=os.path.join(script_destination_folder_on_mesx,new_xml_file_name)
            new_script_file_name = poly_name + ".sh"
            sh_replacements=[("JOBNAME",poly_name),("CONFIGPATH",config_path)]
            new_sh_created=write_config_file(
                template_sh_file,origin_folder,new_script_file_name,sh_replacements)
            self.assertTrue(os.path.exists(new_sh_created))
            cluster_cmds.append(new_script_file_name)

        cmd2=['scp','-r',origin_folder,me_at_remote_URL+':' +commands_folder_on_mesx]
        print(" ".join(cmd2))
        out_string, error_string = process_wrapper.run_and_wait_on_process(cmd2, local_out_dir)

        #can I qsub remotely?
        #ssh "$server" "mkdir -p $destiny"
        #https://stackoverflow.com/questions/26278167/submitting-a-job-to-qsub-remotely
        #https://unix.stackexchange.com/questions/8612/programmatically-creating-a-remote-directory-using-ssh

        cmd3 = ['ssh', me_at_remote_URL, '. ~/.bash_profile;','mkdir ' + specks_output_path_on_mesx]
        print(" ".join(cmd3))
        out_string, error_string = process_wrapper.run_and_wait_on_process(cmd3, local_out_dir)

        for cmd in cluster_cmds:
            cmd4 = ['ssh', me_at_remote_URL, '. ~/.bash_profile;',
                'cd ' + script_destination_folder_on_mesx + ';', 'qsub',cmd]

            print(" ".join(cmd4))
            out_string, error_string = process_wrapper.run_and_wait_on_process(cmd4, local_out_dir)

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
