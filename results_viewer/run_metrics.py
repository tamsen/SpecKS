from results_viewer.curve_fitting import curve_fit_metrics


class run_metrics():
    input_type=''
    sim_name=''
    spc_time=-1
    wgd_time=-1
    lognorm_fit_data=-1
    gaussian_fit_data=-1

    def __init__(self, input_type, sim_name, spc_time, wgd_time,
                 lognorm_fit_data,gaussian_fit_data):
        self.input_type = input_type
        self.sim_name =sim_name
        self.spc_time =spc_time
        self.wgd_time = wgd_time
        self.lognorm_fit_data = lognorm_fit_data
        self.gaussian_fit_data = gaussian_fit_data
                 #mod,cm,num_paralogs,lognorm_RMSE,gaussian_RMSE
    def to_csv_string(self):

        final_data = self.to_data_list()
        return ",".join(final_data)

    def to_data_list(self):
        simple_data = [self.input_type, self.sim_name, self.spc_time, self.wgd_time]#,self.fit_data]
        fit_data_list=self.lognorm_fit_data.to_data_list() + [str(self.gaussian_fit_data.popt)] \
                       +[str(self.gaussian_fit_data.RMSE)]
        final_data = [str(d) for d in simple_data] + [str(d) for d in fit_data_list]
        return final_data

def get_metrics_from_data_line(data_list):
    simple_data=data_list[0:4]
    mode=data_list[4]
    cm=data_list[5]
    np=data_list[6]
    ln_popt= data_list[7]
    ln_RMSE=data_list[8]
    g_popt= data_list[9]
    g_RMSE=data_list[10]
    lognorm_data=curve_fit_metrics(cm,mode,np,ln_popt,ln_RMSE)
    gaussian_data=curve_fit_metrics(cm,mode,np,g_popt,g_RMSE)
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


def get_metric_result_data_headers():

    return ["sim_type","sim_name","spc_time","wgd_time","mode","cm","num_paralogs",
                     "lognorm_popt","lognorm_RMSE","gaussian_popt","gaussian_RMSE"]