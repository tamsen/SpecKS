import math
import os
import unittest
from matplotlib import pyplot as plt
from results_viewer import run_metrics
from results_viewer.batch_analyzer import compute_metrics_for_batch


class BatchAggregator(unittest.TestCase):

    def test_parse_agg_results(self):

        #batch_names = ["sim37_N0p1","sim37_N1","sim37_N5","sim36_N10",
        #               "sim37_N20","sim37_N100","sim35_log"]

        batch_names = ["sim37_N0p1","sim37_N1","sim37_N5","sim36_N10",
                       "sim37_N20","sim35_log"]

        reprocess=False
        #batch_names = ["sim37_N1","sim37_N20","sim37_N100"]
        out_folder = "/home/tamsen/Data/Specks_outout_from_mesx"
        plots_to_make=[([get_spec_time,get_lognorm_RMSE],"spec time (MYA)","lognormRMSE","lognormRMSE.png"),
                       ([get_spec_time,get_gaussian_RMSE], "spec time (MYA)","guassianRMSE","guassianRMSE.png"),
                       ([get_spec_time, get_genes_shed], "spec time (MYA)","# genes shed","genes_shed_vs_spec_time.png"),
                       ([get_wgd_time, get_genes_remaining], "spec time (MYA)", "# genes remaining",
                        "genes_remaining_vs_wgd_time.png"),

                       ([get_wgd_time, get_genes_shed], "wgd time (MYA)","# genes shed","genes_shed_vs_wgd_time.png"),
                       ([get_spec_time, get_mode], "spec time (MYA)", "mode vs spec time", "mode_vs_spec_time.png"),
                       ([get_spec_time, get_max_d], "spec time (MYA)", "max d vs spec time", "maxd_vs_spec_time.png"),
                       ([get_spec_time, get_metric1], "spec time (MYA)", "allo vs auto metric1", "metric1.png"),
                       ([get_spec_time, get_metric2], "spec time (MYA)", "allo vs auto metric2", "metric2.png"),
                       ([get_spec_time, get_metric3], "spec time (MYA)", "allo vs auto metric3 (fit cm-mode)", "metric3.png"),
                       ([get_spec_time, get_metric5_pearsons_r],"spec time (MYA)", "allo vs auto metric5 (fit pr1-pr2)",
                        "metric5.png"),
                       ([get_spec_time, get_metric6_pearsons_pv], "spec time (MYA)", "allo vs auto metric6 (fit pv ratio)",
                        "metric6.png"),
                       ([get_metric7a, get_metric7b], "Ln Pearson corr_coef", "Gaussian pearson corr_coef",
                        "metric7.png"),
                       ([get_wgd_time, get_mode], "wgd time (MYA)", "mode vs wgd time",
                        "mode_vs_wgd_time.png"),
                       ([get_spec_time, get_max], "spec time (MYA)", "max vs spec time",
                        "max_vs_spec_time.png"),
                       ([get_spec_time, get_metric4], "spec time (MYA)",
                        "allo vs auto metric4 (raw cm-mode)", "metric4.png"),
                       ([get_spec_time, get_metric8_pearsons_ratio],"spec time (MYA)",
                       "Pearson corr_coef ratio","metric8.png"),
                        ]
                       #([get_wgd_time, get_mode, get_max_RMSE], "wgd time (MYA)", "mode vs wgd time",
                       # "mode_vs_wgd_time.png"),
                       #([get_spec_time, get_max, get_max_RMSE], "spec time (MYA)", "max vs spec time",
                       # "max_vs_spec_time.png"),
                       #([get_spec_time, get_metric4, get_max_RMSE], "spec time (MYA)",
                       # "allo vs auto metric4 (raw cm-mode)", "metric4.png"),


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
            plot_name = plot_to_make[3]
            plot_data=plot_allo_vs_auto_metrics(metric_data_by_batch,marker_styles_by_batch,
                                                *plot_to_make,out_folder)
            if "metric" in plot_name:
                log_plot_to_make=[d for d in plot_to_make]
                log_plot_to_make[3] =plot_name.replace(".png","_log.png")
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
def sort_batch_into_xs_and_ys_2(metric_data_for_batch,methods):

    allo_xs=[]
    auto_xs=[]
    sim_data=[]

    for sim, run_metrics in metric_data_for_batch.items():

        new_data_point=[]
        for method in methods:
            new_data_point.append(method(run_metrics))

        sim_data.append([sim,*new_data_point])

        if run_metrics.input_type.upper()=="ALLO":
            allo_xs.append(new_data_point)
        else:
            auto_xs.append(new_data_point)

    return allo_xs,auto_xs,sim_data

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

def get_max_RMSE(run_metrics):
    ln_e = get_lognorm_RMSE(run_metrics)
    ga_e = get_gaussian_RMSE(run_metrics)
    return max(ln_e,ga_e)

def get_max_PRSE(run_metrics):
    ln_pearsons_se=float(run_metrics.lognorm_fit_data.get_pearsons_std_error())
    ga_pearsons_se=float(run_metrics.gaussian_fit_data.get_pearsons_std_error())
    return max(ln_pearsons_se,ga_pearsons_se)
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

    #print(run_metrics.sim_name)
    ln_mode=run_metrics.lognorm_fit_data.mode
    ln_cm=run_metrics.lognorm_fit_data.cm
    m3 = ln_cm-ln_mode
    return m3

def get_metric4(run_metrics):

    cm=float(run_metrics.wgd_raw_cm)
    max=float(run_metrics.wgd_raw_x_value_of_ymax)
    return cm-max

def get_metric7a(run_metrics):

    ln_pearsons_cc=float(run_metrics.lognorm_fit_data.get_pearsons_corr_coef())
    return ln_pearsons_cc

def get_metric7b(run_metrics):

    ga_pearsons_cc=float(run_metrics.gaussian_fit_data.get_pearsons_corr_coef())
    return ga_pearsons_cc
def get_metric5_pearsons_r(run_metrics):

    ln_pearsons_cc=float(run_metrics.lognorm_fit_data.get_pearsons_corr_coef())
    ga_pearsons_cc=float(run_metrics.gaussian_fit_data.get_pearsons_corr_coef())
    return ln_pearsons_cc-ga_pearsons_cc

def get_metric8_pearsons_ratio(run_metrics):

    ln_pearsons_cc=float(run_metrics.lognorm_fit_data.get_pearsons_corr_coef())
    ga_pearsons_cc=float(run_metrics.gaussian_fit_data.get_pearsons_corr_coef())
    return math.log( ln_pearsons_cc / ga_pearsons_cc )
def get_metric6_pearsons_pv(run_metrics):

    ln_pearsons_pv=float(run_metrics.lognorm_fit_data.get_pearsons_p_value())
    ga_pearsons_pv=float(run_metrics.gaussian_fit_data.get_pearsons_p_value())

    if (ln_pearsons_pv == 0):
        return math.nan

    log_ga=math.log(ga_pearsons_pv)
    log_pr=math.log(ln_pearsons_pv)
    #ratio=ga_pearsons_pv/ln_pearsons_pv
    return log_ga - log_pr
def get_genes_shed(run_metrics):
    pairs_remaining=run_metrics.lognorm_fit_data.num_paralogs
    num_paralogs_at_WGD=3000
    genes_lost=num_paralogs_at_WGD-pairs_remaining
    y = genes_lost
    return y

def get_genes_remaining(run_metrics):
    return run_metrics.lognorm_fit_data.num_paralogs

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
    sims_without_clear_wgd=[]
    with open(csv_file, "r") as f:

        while True:
            line = f.readline()
            if "sim_type" in line:
                continue
            if len(line)==0:
                break

            data = line.strip().split(",")
            #print(str(data))

            if "NA" in line:
                sims_without_clear_wgd.append(data[1])
                continue

            sim_type=data[0]
            if sim_type=='':
                continue

            metrics_for_run= run_metrics.get_metrics_from_data_line(data)
            run_name=metrics_for_run.sim_name
            data_by_sim_name[run_name]=metrics_for_run

    return data_by_sim_name, sims_without_clear_wgd


def plot_allo_vs_auto_metrics(metric_data_by_batch,marker_styles_by_batch,
                              x_and_y_methods,x_axis_name,y_axis_name,
                              title,out_folder, ylog=False):

    fig, ax = plt.subplots(1, 1, figsize=(8, 10))
    auto_color='gray'
    plot_data={}
    allo_colors = ['cyan','lightblue','blue','darkblue','slateblue','purple','indigo']
    i=0
    for batch_name, metric_data_for_batch in metric_data_by_batch.items():
        marker_style=marker_styles_by_batch[batch_name]

        #[allo_xs,allo_ys],[auto_xs,auto_ys],sim_data =sort_batch_into_xs_and_ys(metric_data_for_batch,get_x_method,get_y_method)
        allo_data, auto_data, sim_data = sort_batch_into_xs_and_ys_2(metric_data_for_batch,
                                                                     x_and_y_methods)
        data_by_label={"allo": allo_data,"auto" :auto_data}
        for label, data in data_by_label.items():

            label_string=batch_name + "({0})".format(label)
            first_data_point=data[0]

            size=30
            if label=="auto":
                color=auto_color
                size=30
            else:
                color=allo_colors[i]

            xs=[d[0] for d in data]
            ys=[d[1] for d in data]
            if len(first_data_point) ==2:
                plt.scatter(xs,ys,
                    c=color,label=label_string,marker=marker_style, s=size)
            else:
                yerr =[d[2] for d in data]
                plt.errorbar(xs,ys, yerr, c=allo_colors[i],label=label,
                     marker=marker_style, linestyle='')

        i=i+1
        plot_data[batch_name] = sim_data

    out_file_name = os.path.join(out_folder, title)
    fig.suptitle(title)

    if ylog:
        ax.set(yscale='log')

    ax.set(xlabel=x_axis_name)
    ax.set(ylabel=y_axis_name)

    #ax.legend()
    #plt.legend(bbox_to_anchor=(1.04, 1), loc="lower center")
    #plt.tight_layout()
    #plt.legend(bbox_to_anchor=(0, -5),loc="lower left",ncol=3)
    plt.legend( bbox_to_anchor=(0,-0.3),loc="lower left", ncol=3)
    fig.set_tight_layout(True)
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

