import os
import unittest
import process_wrapper

class MyTestDownloader(unittest.TestCase):

    def test_ls_results(self):
        me_at_remote_URL = 'tdunn@mesx.sdsu.edu'
        local_output_folder = "/home/tamsen/Data/Specks_outout_from_mesx"
        remote_output_folder = "/usr/scratch2/userdata2/tdunn/SpecKS_Output"
        self.get_subdirectories(local_output_folder, me_at_remote_URL, remote_output_folder)

    def get_subdirectories(self, local_output_folder, me_at_remote_URL, remote_output_folder):
        cmd4 = ['ssh', me_at_remote_URL, '. ~/.bash_profile;',
                'cd ' + remote_output_folder + ';', 'ls']
        print(" ".join(cmd4))
        out_string, error_string = process_wrapper.run_and_wait_on_process(cmd4, local_output_folder)
        return out_string.split()

    def test_download_mesx_results(self):
        batch_folder = "sim33_log"
        #runs = range(1, 7)
        me_at_remote_URL = 'tdunn@mesx.sdsu.edu'
        local_output_folder = "/home/tamsen/Data/Specks_outout_from_mesx"
        remote_output_folder = "/usr/scratch2/userdata2/tdunn/SpecKS_Output"
        local_batch_folder = os.path.join(local_output_folder, batch_folder)
        remote_batch_folder = os.path.join(remote_output_folder, batch_folder)

        if not os.path.exists(local_batch_folder):
            os.makedirs(local_batch_folder)

        run_folders= self.get_subdirectories(local_output_folder, me_at_remote_URL, remote_batch_folder)
        print(run_folders)
        polyploid_folders_by_name={}
        for run_folder in run_folders:
            full_path=os.path.join(remote_batch_folder, run_folder)
            potential_polyploid_folders = self.get_subdirectories(local_output_folder, me_at_remote_URL, full_path)
            for f in potential_polyploid_folders:
                if not "." in f:
                    polyploid_path=os.path.join(full_path, f)
                    polyploid_folders_by_name[f]=polyploid_path

        print(polyploid_folders_by_name)

        for polyploid, remote_path in polyploid_folders_by_name.items():

            #allo_run = "Allo" + str(run)
            #auto_run = "Auto" + str(run)
            local_allo_folder = os.path.join(local_batch_folder, polyploid)
            #local_auto_folder = os.path.join(local_batch_folder, auto_run)

            if not os.path.exists(local_allo_folder):
                os.makedirs(local_allo_folder)
            #if not os.path.exists(local_auto_folder):
            #    os.makedirs(local_auto_folder)

            remote_to_match = remote_path + "/*final*/*.csv"
            cmd2 = ['scp', '-r', me_at_remote_URL + ':' + remote_to_match, local_allo_folder]
            print(" ".join(cmd2))
            #out_string, error_string = process_wrapper.run_and_wait_on_process(cmd2, local_allo_folder)

            #remote_to_match = remote_path + "/*randomized*/*.*"
            #cmd2 = ['scp', '-r', me_at_remote_URL + ':' + remote_to_match, local_allo_folder]
            #print(" ".join(cmd2))
            #out_string, error_string = process_wrapper.run_and_wait_on_process(cmd2, local_batch_folder)


        cmd2 = ['scp', '-r', me_at_remote_URL + ':' + remote_batch_folder + '/specks*/*.xml', '.']
        out_string, error_string = process_wrapper.run_and_wait_on_process(cmd2, local_batch_folder)

        cmd2 = ['scp', '-r', me_at_remote_URL + ':' + remote_batch_folder + '/specks*/*log*', '.']
        out_string, error_string = process_wrapper.run_and_wait_on_process(cmd2, local_batch_folder)


if __name__ == '__main__':
    unittest.main()
