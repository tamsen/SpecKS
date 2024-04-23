import math
import os
import unittest

import numpy as np
from matplotlib import pyplot as plt
import config
import kp_reader
from results_viewer import batch_analyzer, batch_histogrammer, curve_fitting, batch_aggregator


class Test1KP(unittest.TestCase):

    def test_1KP_analysis(self):
        data_directory = "/home/tamsen/Data/1KP_classfier"
        Ks_file_1 = "final_ks_values_CSUV.fa"
        kp_directory = os.path.join(data_directory, "1KP_final_ks_files")
        output_folder = os.path.join(data_directory, "1KP_classifier_out_Apr21")

        lookup_file = "1KP-Sample-List.csv"
        sample_lookup = kp_reader.get_sample_code_lookup(lookup_file)

        if not os.path.exists(output_folder):
            os.makedirs(output_folder)

        reanalyze=True
        bin_size = 0.002#0.02
        kernel_size = 25
        max_Ks = 1.0
        KS_data_files = [Ks_file_1 ] + [f for f in os.listdir(kp_directory) if ".fa" in f]
        max_to_process=50
        i=0
        curated_samples= ['OBUY']#known_auto_examples+known_allo_examples+required
        results_by_file={}

        if not reanalyze:
            KS_data_files=[]

        for i in range(0,len(KS_data_files)):

            Ks_file = KS_data_files[i]
            if i > max_to_process:
                break


            Ks_file_path = os.path.join(kp_directory, Ks_file)
            (species_code, species) = kp_reader.get_species_code(Ks_file, sample_lookup)
            print("classifying " + str(i) + " of " + str(max_to_process))
            ks_values = kp_reader.get_Ks_from_file(Ks_file_path)

            #get hist data
            sample_name=species
            fig, ax = plt.subplots(1, 1, figsize=(10, 10))
            fig.suptitle(sample_name)
            fig.suptitle(sample_name)
            ax.set_title("Allopolyploid\n", fontsize=20)
            ax.set_title("Autopolyploid\n", fontsize=20)

            params=config.PolyploidParams(-1, -1, sample_name)
            hist_data=batch_histogrammer.make_histogram_subplot(ax,  sample_name,ks_values,
                                        bin_size, params,
                                        max_Ks, False,"blue")
            hist_data2=[hist_data[0][1:],hist_data[1][1:]]
            out_file_name = os.path.join(output_folder, "histogram" + "_plot_" +  sample_name+
                                         "_" + species_code + "_"+ str(max_Ks) + ".png")
            plt.savefig(out_file_name)
            plt.close()

            csv_file_name = os.path.join(output_folder, sample_name + ".hist.csv")
            batch_histogrammer.save_hist_to_csv(hist_data2, csv_file_name)

            #analyze hist data
            hist_data3 =  batch_analyzer.read_hist_csv(csv_file_name)
            out_fit_png = os.path.join(output_folder, "fit" + "_plot_" +  sample_name+
                                         "_" + species_code + "_"+ str(max_Ks) + ".png")
            fit_results = batch_analyzer.analyze_histogram(*hist_data3, params.WGD_time_MYA, params.SPC_time_MYA,
                                            kernel_size, False, 'tan', out_fit_png)

            results_by_file[Ks_file]=fit_results
            print(Ks_file + " analyzed")

        out_csv = "1KP_metrics.csv"
        csv_file_name = os.path.join(output_folder, out_csv)
        if reanalyze:
            batch_analyzer.plot_and_save_metrics(results_by_file, csv_file_name)

        analyze_metric_results(csv_file_name,output_folder)

def analyze_metric_results(input_metrics_file, out_folder):

        batch_name = "1KP"
        plots_to_make=[([get_classification, get_log_metric3],
                        "categories", "allo vs auto metric3 (fit cm-mode)",
                        "1KP_allo_vs_auto_classification.png")]
        metric_data_for_batch,sims_without_clear_wgd = batch_aggregator.read_data_csv(input_metrics_file)

        for plot_to_make in plots_to_make:
            x_method=plot_to_make[0][0]
            y_method=plot_to_make[0][1]
            plot_name = plot_to_make[3]
            species_with_WGD_detected_in_1KP=list(metric_data_for_batch.keys())
            plot_data = [[k,x_method(metric_data_for_batch[k]),y_method(metric_data_for_batch[k]) ] for k in
                         species_with_WGD_detected_in_1KP]
            plot_data= plot_data + [[s,"No WGD detected","NA"] for s in sims_without_clear_wgd]

            print(plot_data)
            data_file_name=plot_name.replace(".png",".csv")
            out_file_path=os.path.join(out_folder,data_file_name)
            plot_data_by_batch = {batch_name: plot_data}
            batch_aggregator.save_metrics_to_csv(plot_data_by_batch, out_file_path)
            colors_by_category = {"Low": "gray", "Medium": "cyan", "High": "blue"}
            categories=colors_by_category.keys()
            data=[]
            data_labels=[]
            colors=[]
            for c in categories:
                test_for_data=[d[2] for d in plot_data if d[1]==c]
                if len(test_for_data)> 0:
                    data.append(test_for_data)
                    data_labels.append(c)
                    colors.append(colors_by_category[c])

            ticks=[i+1 for i in range(0,len(data))]
            #[d[2] for d in plot_data if d[1] == "Medium"],
            #https://matplotlib.org/stable/gallery/statistics/violinplot.html
            #https://stackoverflow.com/questions/33864578/matplotlib-making-labels-for-violin-plots
            #colors_by_category = {"Low": "gray", "Medium": "cyan", "High": "blue"}
            #colors=["gray","blue"]

            fig, ax = plt.subplots(1, 1, figsize=(8, 10))
            plots=plt.violinplot(data, ticks,showmeans=True, showextrema=True,
                                 )
            for pc, color in zip(plots['bodies'], colors):
                pc.set_facecolor(color)

            ax.axhline(y=get_low_to_medium_threshold(), color='cyan', linestyle='--', label="lvm disc. criteria"
                                                                         + " ({0})".format(
                round(get_low_to_medium_threshold(), 2)))

            ax.axhline(y=get_medium_to_high_threshold(), color='cyan', linestyle='--', label="lvm disc. criteria"
                                                                         + " ({0})".format(
                round(get_medium_to_high_threshold(), 2)))

            fig.suptitle(plot_to_make[2])
            ax.set(xlabel=plot_to_make[1])
            ax.set(ylabel="log(fit cm-mode)")

            ax.set_xticks(ticks)
            ax.set_xticklabels(data_labels)

            #labels = ['Low', 'High']
            #ax.set_xticks(labels)
            #plt.legend()
            violin_plot_file_path = os.path.join(out_folder, plot_name)
            plt.savefig(violin_plot_file_path)
            plt.close()

def get_low_to_medium_threshold():
    return -4.53

def get_medium_to_high_threshold():
    return -3.12

def get_log_metric3(run_metrics):

    metric_3_result = batch_aggregator.get_metric3(run_metrics)
    if metric_3_result <= 0:
        return math.nan

    return math.log(metric_3_result)
def get_classification(run_metrics):

    metric_3_result_as_log = get_log_metric3(run_metrics)
    low_to_medium_threshold=get_low_to_medium_threshold()
    medium_to_high_threshold=get_medium_to_high_threshold()

    if metric_3_result_as_log < low_to_medium_threshold:
        return "Low"
    if metric_3_result_as_log < medium_to_high_threshold:
        return "Medium"
    if metric_3_result_as_log >= medium_to_high_threshold:
        return "High"
    #else
    return math.nan

if __name__ == '__main__':
    unittest.main()
