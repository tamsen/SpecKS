import os
import unittest

from config import PolyploidParams


class Generate_Config_Files(unittest.TestCase):
    def test_making_configs(self):

        template_dat_file="mesx-template.xml"
        out_dir="/home/tamsen/Data/SpecKS_input"
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)

        decimals_needed=3
        formatter = "{:0" + str(decimals_needed) + "d}"
        spec_times=[80,60,40,20,10]
        wgd_times = [75, 55, 35, 15, 5]
        poly_params_by_name={}

        for i in range(0,len(spec_times)):

            spec_time=spec_times[i]
            wgd_time=wgd_times[i]

            allo_name="Allo"+str(i+1)+"_S"+formatter.format(spec_time)+"W"+formatter.format(wgd_time)
            auto_name="Auto"+str(i+1)+"_S"+formatter.format(spec_time)+"W"+formatter.format(spec_time)

            allo_params = PolyploidParams(spec_time, wgd_time, allo_name)
            auto_params = PolyploidParams(spec_time, spec_time, auto_name)

            poly_params_by_name[allo_name]=allo_params
            poly_params_by_name[auto_name]=auto_params

        for poly_name,poly_params in poly_params_by_name.items():

            new_file_name = poly_name + ".xml"
            new_file_created=write_config_file(
                template_dat_file, out_dir,new_file_name,poly_params)

            self.assertTrue(os.path.exists(new_file_created))

def write_config_file(template_dat_file, out_dir,new_file_name,poly_params):

    lines_to_write = []
    new_config_file = os.path.join(out_dir, new_file_name)

    with open(template_dat_file, 'r') as f:

        while (True):

            line = f.readline()
            new_line = line

            if "POLYPLOID_SECTION" in line:
                new_line = line.replace("POLYPLOID_SECTION", poly_params.to_xml())


            lines_to_write.append(new_line)

            if not line:
                break

    with open(new_config_file, 'w') as f:

        for line in lines_to_write:
            f.writelines(line)

    return new_config_file



if __name__ == '__main__':
    unittest.main()
