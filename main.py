import os
import allosim
import autosim
import config
from datetime import datetime

from polyploid_setup import make_polyploids


def run_sim():

    conf = setup()

    # Time since WGD: 5,10, 15,50,100,200 MYA. Total tree length 500 MY. Make allo and autopoly examples.
    list_of_polyploids = make_polyploids(conf)

    for polyploid in list_of_polyploids[1:2]:

        if polyploid.is_allo():
            allosim.run_allosim(polyploid)
        else:
            autosim.run_autosim(polyploid)

    print("\n\nSpecKS complete")

def setup():

    now = datetime.now()
    date_time = now.strftime("m%md%dy%Y_h%Hm%Ms%S")
    conf = config.SpecKS_config()
    conf.output_folder = conf.output_folder_root + "_" + date_time

    print(conf.output_folder)
    if not os.path.exists(conf.output_folder):
        os.makedirs(conf.output_folder)

    print("Current environment:")
    print(str(os.environ))

    cwd=os.getcwd()
    print("Current Working Directory:\t" + cwd)

    return conf


if __name__ == '__main__':
    run_sim()

# See PyCharm help at https://www.jetbrains.com/help/pycharm/
