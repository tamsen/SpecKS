import os
import unittest

from matplotlib import pyplot as plt


class MyAggparser(unittest.TestCase):

    def test_parse_agg_results(self):

        out_folder= "/home/tamsen/Data/Specks_outout_from_mesx/"
        csv_file="/home/tamsen/Data/Specks_outout_from_mesx/joint_data.csv"
        data=read_data_csv(csv_file)

        #get_allo_vs_auto_data(out_folder, xs, ys, title):
        xs_by_sim=[]
        ys_by_sim=[]
        labels_by_sim=[]
        colors_by_sim=[]
        type_by_sim=[]
        for sim, sim_result in data.items():
            for polyploid_type, polyploid_results in sim_result.items():
                    xs=[res.spc_time  for res in polyploid_results]
                    xs_by_sim.append(xs)
                    ys = [abs(res.cm-res.mode) for res in polyploid_results]
                    ys_by_sim.append(ys)
                    labels_by_sim.append(sim)
                    type_by_sim.append(polyploid_type)
                    if polyploid_type== "Allo":
                        colors_by_sim.append('r')
                    else:
                        colors_by_sim.append('b')
        print(data)
        plot_allo_vs_auto_metrics(out_folder,xs_by_sim,ys_by_sim, labels_by_sim, colors_by_sim,type_by_sim,
                                  "Allopolyploids show more skew Ks histograms than Autopolyploids ")
        self.assertEqual(True, False)  # add assertion here


def read_data_csv(csv_file):

    data_by_sim_name_by_type={}
    with open(csv_file, "r") as f:

        while True:
            line = f.readline()
            if "input_type" in line:
                continue
            if len(line)==0:
                break
            data = line.strip().split(",")
            sim_name=data[0]
            if sim_name=='':
                continue

            if sim_name not in data_by_sim_name_by_type:
                data_by_sim_name_by_type[sim_name]={}

            polyploid_type_string=data[1]
            polyploid_type="Auto"
            if "Allo" in polyploid_type_string:
                polyploid_type="Allo"

            if polyploid_type not in data_by_sim_name_by_type[sim_name]:
                data_by_sim_name_by_type[sim_name][polyploid_type]=[]

            sim_metric_result= metric_result(data)
            data_by_sim_name_by_type[sim_name][polyploid_type].append(sim_metric_result)

    return data_by_sim_name_by_type


def plot_allo_vs_auto_metrics(out_folder, list_of_xs, list_of_ys, list_of_sims,colors_by_sim,type_by_sim, title):

    marker_styles={"lognorm":"x","N=0p1":"+","N=1p0" : "o"}
    fig, ax = plt.subplots(1, 1, figsize=(10, 10))

    for i in range(0,len(list_of_xs)):
        label= list_of_sims[i] + "(" + type_by_sim[i] + ")"
        plt.scatter(list_of_xs[i], list_of_ys[i], label=label,c=colors_by_sim[i],
                    marker=marker_styles[list_of_sims[i]] )
    out_file_name = os.path.join(out_folder, title+ ".png")
    fig.suptitle(title)
    ax.set(xlabel="Spec time, MYA")
    ax.set(ylabel="Metric (mode-cm)")

    ax.legend()
    plt.savefig(out_file_name)
    plt.close()
if __name__ == '__main__':
    unittest.main()

class metric_result():
    input_type=''
    sim_name=''
    spc_time=-1
    wgd_time=-1
    mode=-1
    cm=-1
    num_paralogs=-1
    popt=[]
    norm_fit_result=[]
    lognorm_fit_result=[]
    def __init__(self,data_list):
        self.input_type = data_list[0]
        self.sim_name =data_list[1]
        self.spc_time =int(data_list[2])
        self.wgd_time = int(data_list[3])
        self.mode = float(data_list[4])
        self.cm = float(data_list[5])
        self.num_paralogs = int(data_list[6])
        self.popt = [float(d) for d in data_list[7:7+4]]
        self.norm_fit_result=data_list[11]
        self.lognorm_fit_result=data_list[12]
    def to_csv_string(self):

        final_data = self.to_data_list()
        return ",".join(final_data)

    def to_data_list(self):
        simple_data = [self.input_type, self.sim_name, self.spc_time, self.wgd_time,
                       self.mode, self.cm, self.num_paralogs]
        p_opt_data = [p for p in self.popt]
        ks_data =[str(self.norm_fit_result), str(self.norm_fit_result)]
        final_data = simple_data + p_opt_data + ks_data
        return final_data

def get_metric_result_data_headers():
         return ["sim_name","spc_time","wgd_time","mode","cm","num_paralogs",
                 "popt1","popt2","popt3","popt4","norm_fit_result","lognorm_fit_result"]