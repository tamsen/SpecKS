import os
import sys

import allosim
import autosim
import config
from datetime import datetime

from polyploid_setup import make_polyploids


def run_sim():

    conf = setup(sys.argv)

    # Time since WGD: 5,10, 15,50,100,200 MYA. Total tree length 500 MY. Make allo and autopoly examples.
    list_of_polyploids = make_polyploids(conf)

    for polyploid in list_of_polyploids[0:2]:

        if polyploid.is_allo():
            allosim.run_allosim(polyploid)
        else:
            autosim.run_autosim(polyploid)

    print("\n\nSpecKS complete")

def setup(arguments):

    config_file=arguments[1]
    now = datetime.now()
    date_time = now.strftime("m%md%dy%Y_h%Hm%Ms%S")
    conf = config.SpecKS_config(config_file)
    conf.output_folder = conf.output_folder_root + "_" + date_time

    print('Command Arguments Given: %s' % arguments)
    print('Conifg file: %s' % config_file)
    print("Current environment: %s" + str(os.environ))
    print("Current Working Directory:\t" + os.getcwd())

    print("Output filder:\t" + conf.output_folder)
    if not os.path.exists(conf.output_folder):
        os.makedirs(conf.output_folder)

    return conf


if __name__ == '__main__':
    run_sim()

# See PyCharm help at https://www.jetbrains.com/help/pycharm/
