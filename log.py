import os
from datetime import datetime
import version

log_file_path = ""

def write_start_to_log( outputDir, log_file_name, version_info ):

    global log_file_path
    log_file_path = os.path.join(outputDir,log_file_name)
    write_to_log("\n")
    write_to_log("*******************************************")
    write_to_log("************** SpecKS initiated ***********")

    version_line = version_info.to_string()
    write_to_log(version_line)
    write_to_log("*******************************************")

def write_end_to_log( ):

    write_to_log("*********** SpecKS completed **************")
    write_to_log("*******************************************")
    write_to_log("\n")

def write_error_to_log(msg):

    write_to_log("*******************Error*******************")
    write_to_log("Error: " + msg)


def write_warning_to_log(msg):
    write_to_log("*******Warning******")
    write_to_log("Warning*: " + msg)
    write_to_log(msg)

def write_to_log(msg):

        now = datetime.now()
        day = now.strftime("%d/%m/%Y")
        time = now.strftime("%H:%M:%S")
        time_stamp_string = ",".join([day,time])
        log_line = time_stamp_string + ":\t" + msg + "\n"
        print(log_line)

        if len(log_file_path) == 0:
            return

        with open(log_file_path, 'a') as f:
           f.write(log_line)