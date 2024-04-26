from results_viewer import curve_fitting
from results_viewer.curve_fitting import curve_fit_metrics


class run_metrics():
    input_type=''
    sim_name=''
    spc_time=-1
    wgd_time=-1
    wgd_raw_cm=-1
    wgd_raw_x_value_of_ymax=-1
    wgd_fit_maxima=-1
    wgd_fit_maxima_d=-1
    lognorm_fit_data=-1
    gaussian_fit_data=-1

    def __init__(self, input_type, sim_name, spc_time, wgd_time,
                 raw_cm, raw_x_value_of_ymax,
                 wgd_fit_maxima, wgd_fit_maxima_d, lognorm_fit_data, gaussian_fit_data):
        self.input_type = input_type
        self.sim_name =sim_name
        self.spc_time =spc_time
        self.wgd_time = wgd_time
        self.wgd_raw_cm = raw_cm
        self.wgd_raw_x_value_of_ymax =raw_x_value_of_ymax
        self.wgd_fit_maxima = wgd_fit_maxima
        self.wgd_fit_maxima_d = wgd_fit_maxima_d
        self.lognorm_fit_data = lognorm_fit_data
        self.gaussian_fit_data = gaussian_fit_data
                 #mod,cm,num_paralogs,lognorm_RMSE,gaussian_RMSE
    def to_csv_string(self):

        final_data = self.to_data_list()
        return ",".join(final_data)

    def to_data_list(self):
        simple_data = [self.input_type, self.sim_name, self.spc_time,
                       self.wgd_time, self.wgd_raw_cm, self.wgd_raw_x_value_of_ymax,
                       self.wgd_fit_maxima, self.wgd_fit_maxima_d]

        if self.lognorm_fit_data and self.gaussian_fit_data:
            fit_data_list=self.lognorm_fit_data.to_data_list() + self.gaussian_fit_data.to_data_list()
        else:
            fit_data_list =  ["NA","NA", "NA","NA", "NA","NA","NA","NA","NA","NA","NA","NA","NA"]
        final_data = [str(d) for d in simple_data] + [str(d) for d in fit_data_list]
        return final_data

def get_metric_result_data_headers():

    ln_headers= ["ln_" + h for h in curve_fitting.col_headers_for_curve_fit_metrics()]
    ga_headers= ["ga_" + h for h in curve_fitting.col_headers_for_curve_fit_metrics()]

    return ["sim_type","sim_name","spc_time",
            "wgd_time","raw_cm","raw_x_value_of_ymax",
            "wgd_fit_maxima","wgd_fit_maxima_d",
            *ln_headers,*ga_headers]
def get_metrics_from_data_line(data_list):

    simple_data=data_list[0:8]
    ln_data=data_list[8:15]
    ga_data=data_list[15:22]
    lognorm_data=curve_fit_metrics(*ln_data)
    gaussian_data=curve_fit_metrics(*ga_data)
    metrics=run_metrics(*simple_data,lognorm_data,gaussian_data)
    return metrics
def plot_and_save_metrics(metrics_by_result_names, out_file_name):

    allo_results={}
    auto_results={}
    for sim_name, metric in metrics_by_result_names.items():
        if "Allo" in sim_name:
            allo_results[sim_name]=metric
        else:
            auto_results[sim_name] = metric

    save_metrics_to_csv(allo_results, auto_results, out_file_name)


def save_metrics_to_csv(allo_results, auto_results, out_file_name):

    with open(out_file_name, 'w') as f:
        metric_result_data_headers= get_metric_result_data_headers()
        f.writelines(",".join(metric_result_data_headers) +"\n")

        for sim_name, metric in allo_results.items():
            f.writelines(metric.to_csv_string() + "\n")

        for sim_name, metric in auto_results.items():
            f.writelines(metric.to_csv_string() + "\n")

def read_run_metrics_from_csv(csv_file_name):

    with open(csv_file_name, 'r') as f:

        while True:
            line = f.readline()
            if "sim_type" in line:
                continue
            if len(line) == 0:
                break
            data = line.strip().split(",")
            sim_name = data[0]
            if sim_name == '':
                continue

