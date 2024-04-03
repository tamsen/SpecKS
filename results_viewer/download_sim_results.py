import os
import unittest
import process_wrapper

class MyTestDownloader(unittest.TestCase):

    def test_download_mesx_results(self):
        batch_folder = "sim26_log"
        runs = range(1, 7)
        me_at_remote_URL = 'tdunn@mesx.sdsu.edu'
        local_output_folder = "/home/tamsen/Data/Specks_outout_from_mesx"
        remote_output_folder = "/usr/scratch2/userdata2/tdunn/SpecKS_Output"
        local_batch_folder = os.path.join(local_output_folder, batch_folder)
        remote_batch_folder = os.path.join(remote_output_folder, batch_folder)

        if not os.path.exists(local_batch_folder):
            os.makedirs(local_batch_folder)

        for run in runs:
            allo_run = "Allo" + str(run)
            auto_run = "Auto" + str(run)
            local_allo_folder = os.path.join(local_batch_folder, allo_run)
            local_auto_folder = os.path.join(local_batch_folder, auto_run)

            if not os.path.exists(local_allo_folder):
                os.makedirs(local_allo_folder)
            if not os.path.exists(local_auto_folder):
                os.makedirs(local_auto_folder)

            remote_allo_folder = os.path.join(remote_batch_folder, "specks_" + allo_run.upper() + "*")
            remote_to_match = remote_allo_folder + "/A*/*final*/*.csv"
            cmd2 = ['scp', '-r', me_at_remote_URL + ':' + remote_to_match, local_allo_folder]
            print(" ".join(cmd2))
            out_string, error_string = process_wrapper.run_and_wait_on_process(cmd2, local_batch_folder)

            remote_auto_folder = os.path.join(remote_batch_folder, "specks_" + auto_run.upper() + "*")
            remote_to_match = remote_auto_folder + "/A*/*final*/*.csv"
            cmd2 = ['scp', '-r', me_at_remote_URL + ':' + remote_to_match, local_auto_folder]
            print(" ".join(cmd2))
            out_string, error_string = process_wrapper.run_and_wait_on_process(cmd2, local_batch_folder)

        cmd2 = ['scp', '-r', me_at_remote_URL + ':' + remote_batch_folder + '/specks*/*.xml', '.']
        out_string, error_string = process_wrapper.run_and_wait_on_process(cmd2, local_batch_folder)

        cmd2 = ['scp', '-r', me_at_remote_URL + ':' + remote_batch_folder + '/specks*/*log*', '.']
        out_string, error_string = process_wrapper.run_and_wait_on_process(cmd2, local_batch_folder)


if __name__ == '__main__':
    unittest.main()
