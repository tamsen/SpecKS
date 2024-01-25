import os
import subprocess
def run_and_wait_on_process(cmd, folder):

    program=cmd[0]
    result_CompletedProcess = subprocess.run(cmd, capture_output=True, cwd=folder)
    error_string=result_CompletedProcess.stderr.decode()
    out_string=result_CompletedProcess.stdout.decode()

    #https://stackoverflow.com/questions/287871/how-do-i-print-colored-text-to-the-terminal
    #colored_error_string = '\x1b[6;30;42m' + "ERROR:  " + error_string + '\x1b[0m'
    colored_error_string = '\033[93m' + "ERROR:  " + error_string + '\x1b[0m'

    if len(error_string)> 0:
        print(colored_error_string )

    #print(program+" stderr:\t" + error_string)
    #print(program+" stdout:\t" + out_string)

    with open(os.path.join(folder, program+ "_stderr.txt"), 'w') as f:
        f.writelines(error_string)
    with open(os.path.join(folder, program+"_stdout.txt"), 'w') as f:
        f.writelines(out_string)