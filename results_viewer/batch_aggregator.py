import os
import unittest
from matplotlib import pyplot as plt
from results_viewer import run_metrics
from results_viewer.batch_analyzer import compute_metrics_for_batch


class BatchAggregator(unittest.TestCase):

    def test_parse_agg_results(self):

        batch_names = ["sim37_N0p1","sim37_N1","sim37_N5","sim36_N10",
                       "sim37_N20","sim37_N100","sim35_log"]

        out_folder = "/home/tamsen/Data/Specks_outout_from_mesx"
        plots_to_make=[(get_spec_time,get_lognorm_RMSE,"spec time (MYA)","lognormRMSE","lognormRMSE.png"),
                       (get_spec_time,get_gaussian_RMSE, "spec time (MYA)","guassianRMSE","guassianRMSE.png"),
                       (get_spec_time, get_genes_shed, "spec time (MYA)","# genes shed","genes_shed_vs_spec_time.png"),
                       (get_wgd_time, get_genes_shed, "wgd time (MYA)","# genes shed","genes_shed_vs_wgd_time.png"),
                       (get_spec_time, get_mode, "spec time (MYA)", "mode vs spec time", "mode_vs_spec_time.png"),
                       (get_wgd_time, get_mode, "wgd time (MYA)", "mode vs wgd time", "mode_vs_wgd_time.png"),
                       (get_spec_time, get_max, "spec time (MYA)", "max vs spec time", "max_vs_spec_time.png"),
                       (get_spec_time, get_max_d, "spec time (MYA)", "max d vs spec time", "maxd_vs_spec_time.png"),
                       (get_spec_time, get_metric4, "spec time (MYA)", "allo vs auto metric4 (raw cm-mode)", "metric4.png"),
                       (get_spec_time, get_metric1, "spec time (MYA)", "allo vs auto metric1", "metric1.png"),
                       (get_spec_time, get_metric2, "spec time (MYA)", "allo vs auto metric2", "metric2.png"),
                       (get_spec_time, get_metric3, "spec time (MYA)", "allo vs auto metric3 (fit cm-mode)", "metric3.png"),
                       ]

        reprocess=False
        marker_styles_for_batches = [".", "+", "*", ">","<", "^", "*", ">"]

        if reprocess:
            for batch_name in batch_names:
                reprocess_batch(batch_name, out_folder)

        metric_data_by_batch = self.read_in_metric_data_fo_all_batches(batch_names, out_folder)

        i=0; marker_styles_by_batch={}
        for batch_name in batch_names:
            marker_styles_by_batch[batch_name]=marker_styles_for_batches[i]
            i=i+1

        for plot_to_make in plots_to_make:
            plot_name = plot_to_make[4]
            plot_data=plot_allo_vs_auto_metrics(metric_data_by_batch,marker_styles_by_batch, *plot_to_make,out_folder)
            if "metric" in plot_name:
                log_plot_to_make=[d for d in plot_to_make]
                log_plot_to_make[4] =plot_name.replace(".png","_log.png")
                plot_allo_vs_auto_metrics(metric_data_by_batch, marker_styles_by_batch,
                                          *log_plot_to_make, out_folder, ylog=True)
            data_file_name=plot_name.replace(".png",".csv")
            out_file_path=os.path.join(out_folder,data_file_name)
            save_metrics_to_csv(plot_data, out_file_path)

    def read_in_metric_data_fo_all_batches(self, batch_names, out_folder):
        metric_data_by_batch = {}
        for batch_name in batch_names:
            batch_folder = os.path.join(out_folder, batch_name, "analysis")
            files = os.listdir(batch_folder)
            for file in files:
                print("file:\t" + file)
                if ("metrics.csv" in file):# and ("polyploid" in file) and ("batch_processed" in file):
                    full_file_path = os.path.join(batch_folder, file)
                    print("full_file_path:\t" + full_file_path)
                    metric_data_for_batch = read_data_csv(full_file_path)
                    metric_data_by_batch[batch_name] = metric_data_for_batch
        return metric_data_by_batch


def reprocess_batch(batch_name,output_folder):
    compute_metrics_for_batch(batch_name,output_folder)
def sort_batch_into_xs_and_ys(metric_data_for_batch,get_x,get_y):

    allo_xs=[]
    allo_ys=[]
    auto_xs=[]
    auto_ys=[]
    sim_data=[]

    for sim, run_metrics in metric_data_for_batch.items():

        x = get_x(run_metrics)
        y = get_y(run_metrics)
        sim_data.append([sim,x,y])

        if run_metrics.input_type.upper()=="ALLO":
            allo_xs.append(x)
            allo_ys.append(y)
        else:
            auto_xs.append(x)
            auto_ys.append(y)

    return [allo_xs,allo_ys],[auto_xs,auto_ys],sim_data


def get_lognorm_RMSE(run_metrics):
    y = round(float(run_metrics.lognorm_fit_data.RMSE),3)
    return y

def get_gaussian_RMSE(run_metrics):
    y = round(float(run_metrics.gaussian_fit_data.RMSE),3)
    return y



def get_metric1(run_metrics):

    ln_rmse=get_lognorm_RMSE(run_metrics)
    g_rmse=get_gaussian_RMSE(run_metrics)

    if g_rmse ==0:
        return 1.0
    m1 = (g_rmse - ln_rmse)/g_rmse
    return m1

def get_metric2(run_metrics):

    ln_rmse=get_lognorm_RMSE(run_metrics)
    g_rmse=get_gaussian_RMSE(run_metrics)

    m1 = (g_rmse - ln_rmse)/1.0
    return m1

def get_metric3(run_metrics):

    ln_mode=run_metrics.lognorm_fit_data.mode
    ln_cm=run_metrics.lognorm_fit_data.cm
    m3 = ln_cm-ln_mode
    return m3

def get_metric4(run_metrics):

    cm=float(run_metrics.wgd_raw_cm)
    max=float(run_metrics.wgd_raw_x_value_of_ymax)
    return cm-max
def get_genes_shed(run_metrics):
    pairs_remaining=run_metrics.lognorm_fit_data.num_paralogs
    num_paralogs_at_WGD=3000
    genes_lost=num_paralogs_at_WGD-pairs_remaining
    y = genes_lost
    return y

def get_spec_time(run_metrics):
    return int(run_metrics.spc_time)

def get_wgd_time(run_metrics):
    return int(run_metrics.wgd_time)

def get_mode(run_metrics):
    y = round(float(run_metrics.lognorm_fit_data.mode),3)
    return y

def get_max(run_metrics):
    y = round(float(run_metrics.wgd_fit_maxima), 3)
    return y

def get_max_d(run_metrics):
    y = round(float(run_metrics.wgd_fit_maxima_d), 3)
    return y
def read_data_csv(csv_file):

    data_by_sim_name={}
    with open(csv_file, "r") as f:

        while True:
            line = f.readline()
            if "sim_type" in line:
                continue
            if len(line)==0:
                break
            if "NA" in line:
                continue
            data = line.strip().split(",")
            print(str(data))
            sim_type=data[0]
            if sim_type=='':
                continue

            metrics_for_run= run_metrics.get_metrics_from_data_line(data)
            run_name=metrics_for_run.sim_name
            data_by_sim_name[run_name]=metrics_for_run

    return data_by_sim_name


def plot_allo_vs_auto_metrics(metric_data_by_batch,marker_styles_by_batch,
                              get_x_method,get_y_method,x_axis_name,y_axis_name,
                              title,out_folder, ylog=False):

    fig, ax = plt.subplots(1, 1, figsize=(5, 5))
    auto_color='gray'
    plot_data={}
    allo_colors = ['cyan','lightblue','blue','darkblue','slateblue','purple','indigo']
    i=0
    for batch_name, metric_data_for_batch in metric_data_by_batch.items():
        marker_style=marker_styles_by_batch[batch_name]
        [allo_xs,allo_ys],[auto_xs,auto_ys],sim_data =sort_batch_into_xs_and_ys(metric_data_for_batch,get_x_method,get_y_method)

        label=batch_name + "(allo)"
        plt.scatter(allo_xs,allo_ys,c=allo_colors[i],label=label,marker=marker_style)
        label=batch_name + "(auto)"
        plt.scatter(auto_xs,auto_ys,c=auto_color,label=label,marker=marker_style)

        i=i+1
        plot_data[batch_name] = sim_data

    out_file_name = os.path.join(out_folder, title)
    fig.suptitle(title)

    if ylog:
        ax.set(yscale='log')

    ax.set(xlabel=x_axis_name)
    ax.set(ylabel=y_axis_name)

    ax.legend()
    plt.savefig(out_file_name)
    plt.close()
    return plot_data
def save_metrics_to_csv(plot_data, out_file_name):

    with open(out_file_name, 'w') as f:
        data_headers= ['batch','sim','x','y']
        f.writelines(",".join(data_headers) +"\n")

        for batch_name, plot_data_for_sims in plot_data.items():

            for plot_data in plot_data_for_sims:
                data_list=[batch_name,*plot_data]
                data_list_string=[str(d) for d in data_list]
                f.writelines(",".join(data_list_string) +"\n")


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
        #self.popt = [float(d) for d in data_list[7:7+4]]
        #self.norm_fit_result=data_list[11]
        #self.lognorm_fit_result=data_list[12]
    def to_csv_string(self):

        final_data = self.to_data_list()
        return ",".join(final_data)

    def to_data_list(self):
        simple_data = [self.input_type, self.sim_name, self.spc_time, self.wgd_time,
                       self.mode, self.cm, self.num_paralogs]
        #p_opt_data = [p for p in self.popt]
        #ks_data =[str(self.norm_fit_result), str(self.norm_fit_result)]
        final_data = simple_data #+ p_opt_data + ks_data
        return final_data

def get_metric_result_data_headers():
         return ["sim_name","spc_time","wgd_time","mode","cm","num_paralogs",
                 "popt1","popt2","popt3","popt4","norm_fit_result","lognorm_fit_result"]
