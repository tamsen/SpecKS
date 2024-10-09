import os
import shutil
import sys
import polyploid_sim
import config
from datetime import datetime

import log
import version
from polyploid_setup import make_polyploids


def run_sim():

    conf = setup(sys.argv)
    if not conf:
        return

    #start the log
    log.write_start_to_log(conf.output_folder,conf.log_file_name, conf.version_info)
    log.write_to_log('Command Arguments Given: %s' % sys.argv)


    if not conf.evolver_random_seed: #if not set by the user...
        log.write_to_log('evolver_random_seed not set by user')
        conf.evolver_random_seed = 137
        log.write_to_log('evolver_random_seed set to ' + str(conf.evolver_random_seed)
                         + " by default.")
    else:
        log.write_to_log('evolver_random_seed set to ' + str(conf.evolver_random_seed)
                         + " by user.")

    if not conf.specks_random_seed: #if not set by the user...
        log.write_to_log('specks_random_seed not set by user')
        conf.evolver_random_seed = 84
        log.write_to_log('specks_random_seed set to ' + str(conf.specks_random_seed)
                         + " by default.")
    else:
        log.write_to_log('specks_random_seed set to ' + str(conf.specks_random_seed)
                         + " by user.")

    #run the sim
    list_of_polyploids = make_polyploids(conf)
    for polyploid in list_of_polyploids:
        polyploid_sim.run_sim(polyploid)

    log.write_end_to_log()

def setup(arguments):

    print('Command Arguments Given: %s' % arguments)
    if len(arguments) < 2:
        print('Please give an input file path.')
        return False

    config_file=arguments[1]
    now = datetime.now()
    date_time = now.strftime("m%md%dy%Y_h%Hm%Ms%S")
    conf = config.SpecKS_config(config_file)
    conf.output_folder = conf.output_folder_root + "_" + date_time
    conf.log_file_name = date_time + "_" + conf.log_file_name
    conf.version_info = version.version_info()
    cwd=os.getcwd()

    print('Config file: %s' % config_file)
    print("Current environment: %s" + str(os.environ))
    print("Current Working Directory:\t" + cwd)
    if conf.output_folder[0:2]== "./":
        conf.output_folder = os.path.join(os.getcwd(),conf.output_folder.replace("./",""))

    config_file_used=os.path.basename(config_file).replace(".xml",".used.xml")
    print("Output folder:\t" + conf.output_folder)
    if not os.path.exists(conf.output_folder):
        os.makedirs(conf.output_folder)

    #move a copy of the config file into the output folder so we remember what was run
    dst = os.path.join(conf.output_folder,config_file_used)
    shutil.copyfile(config_file, dst)

    return conf


if __name__ == '__main__':
    run_sim()
