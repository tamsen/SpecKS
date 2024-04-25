import os
import time
import unittest
import process_wrapper

class MyTestDownloader(unittest.TestCase):



    def get_subdirectories(self, local_output_folder, me_at_remote_URL, remote_output_folder, argument=False):
        cmd4 = ['ssh', me_at_remote_URL, '. ~/.bash_profile;',
                'cd ' + remote_output_folder + ';', 'ls']
        if argument:
            cmd4.append(argument)
        print(" ".join(cmd4))
        #out_string, error_string = process_wrapper.run_and_wait_on_process(cmd4, local_output_folder)
        out_string, error_string = process_wrapper.run_and_wait_with_retry(cmd4, local_output_folder,
                                        "Connection reset by peer", 3,5)
        return out_string

    def test_download_mesx_results(self):

        batch_folder = "sim39_0p1" #"sim37_N20" #sim37_N0p1,sim37_N5
        me_at_remote_URL = 'tdunn@mesx.sdsu.edu'
        local_output_folder = "/home/tamsen/Data/Specks_outout_from_mesx"
        remote_output_folder = "/usr/scratch2/userdata2/tdunn/SpecKS_Output"
        local_batch_folder = os.path.join(local_output_folder, batch_folder)
        remote_batch_folder = os.path.join(remote_output_folder, batch_folder)

        if not os.path.exists(local_batch_folder):
            os.makedirs(local_batch_folder)

        run_folders= self.get_subdirectories(local_output_folder, me_at_remote_URL,
                                             remote_batch_folder, "sp*")
        run_folder_by_polyploid_name = self.get_run_folders_by_polyploid_name(run_folders)


        polyploid_folders_by_name={}
        for polyploid_name, run_folder in run_folder_by_polyploid_name.items():
            full_path=os.path.join(remote_batch_folder, run_folder,polyploid_name)
            polyploid_folders_by_name[polyploid_name]=full_path

        print(polyploid_folders_by_name)

        scp_commands=[]
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
            out_string, error_string = process_wrapper.run_and_wait_with_retry(cmd2, local_allo_folder,
                                            "Connection reset by peer", 2, 5)
            time.sleep(5)
            #scp_commands.append(polyploid)
            scp_commands.append(" ".join(cmd2))
            #remote_to_match = remote_path + "/*randomized*/*.*"
            #cmd2 = ['scp', '-r', me_at_remote_URL + ':' + remote_to_match, local_allo_folder]
            #print(" ".join(cmd2))
            #out_string, error_string = process_wrapper.run_and_wait_on_process(cmd2, local_batch_folder)

        time.sleep(5)
        cmd2 = ['scp', '-r', me_at_remote_URL + ':' + remote_batch_folder + '/specks*/*.xml', '.']
        print(" ".join(cmd2))
        out_string, error_string = process_wrapper.run_and_wait_on_process(cmd2, local_batch_folder)

        time.sleep(5)
        cmd2 = ['scp', '-r', me_at_remote_URL + ':' + remote_batch_folder + '/specks*/*log*', '.']
        print(" ".join(cmd2))
        out_string, error_string = process_wrapper.run_and_wait_on_process(cmd2, local_batch_folder)

        print(scp_commands)

    def get_run_folders_by_polyploid_name(self, run_folders):
        directory_groups = run_folders.split("\n\n")
        #print(run_folders)
        #print(directory_groups)
        run_folder_by_polyploid_name = {}
        for group in directory_groups:
            paths = group.split("\n")
            directory_name = paths[0].replace(":", "")
            poyploid_names = []
            for path in paths[1:len(paths)]:

                if path == "":
                    continue

                if not "." in path:
                    poyploid_names.append(path)

            for polyploid in poyploid_names:
                run_folder_by_polyploid_name[polyploid] = directory_name
        return run_folder_by_polyploid_name


if __name__ == '__main__':
    unittest.main()
