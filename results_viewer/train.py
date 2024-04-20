import math
import os
import unittest
import numpy as np
from sklearn import metrics
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import confusion_matrix
from sklearn.model_selection import train_test_split
import matplotlib.pyplot as plt

from results_viewer.allo_auto_predictor import read_xs_ys_csv, categorize_sim


#https://stackoverflow.com/questions/25009284/how-to-plot-roc-curve-in-python
# https://stackoverflow.com/questions/28716241/controlling-the-threshold-in-logistic-regression-in-scikit-learn

class TestLinearRegression(unittest.TestCase):

    def test_linear_regression_on_allo_vs_auto_data(self):

        out_folder = "/home/tamsen/Data/Specks_outout_from_mesx"
        file_basename= "allopolyploid_index_truth_vs_predictions.csv"
        csv_file=os.path.join(out_folder,file_basename)
        specs, data, target = load_allo_vs_auto_data(csv_file)
        colors = ["red" if m == 0 else "gray" for m in target]
        linear_regression_threshold_color = 'black'
        likely_range_for_threshold = np.arange(0.8, 2, 0.001)
        ROC_plot_base_name = "allo_vs_auto_basic_ROC2"
        threshold_plot_title = "threshold_on_allo_vs_auto_data"
        threshold_plot_labels_by_case={0:"Auto",1:"Allo"}
        case0_color="red"
        self.make_both_ROC_plots(ROC_plot_base_name, colors, data, likely_range_for_threshold,
                                 linear_regression_threshold_color, case0_color,
                                 out_folder, specs, target,
                                 threshold_plot_title, threshold_plot_labels_by_case)

    def test_linear_regression_on_high_vs_low_N_data(self):

        out_folder = "/home/tamsen/Data/Specks_outout_from_mesx"
        sims,specs,true_category,data=get_highN_vs_lowN_truth()

        target1 = [0 if m == "Low" else 1 for m in true_category]
        colors1 = ["gray" if m == 0  else "cyan" for m in target1]
        linear_regression_threshold_color1 = 'cyan'
        likely_range_for_threshold = np.arange(-6, -4, 0.001)
        ROC_plot_base_name="low_vs_high&medium_N_basic_ROC"
        threshold_plot_title="threshold_on_low_vs_high&medium_N_data"
        threshold_plot_labels_by_case={0:"Low",1:"Medium&High"}
        case0_color="gray"
        linear_regression_threshold1=make_both_ROC_plots(ROC_plot_base_name,
                                 colors1, data, likely_range_for_threshold,
                                 linear_regression_threshold_color1, case0_color,
                                 out_folder, specs, target1,
                                 threshold_plot_title, threshold_plot_labels_by_case)

        target2 = [1 if m == "High" else 0 for m in true_category]
        colors2 = ["cyan" if m == 0 else "blue" for m in target2]
        linear_regression_threshold_color2 = 'blue'
        likely_range_for_threshold = np.arange(-6, -1, 0.001)
        ROC_plot_base_name = "low&medium_vs_high_N_basic_ROC"
        threshold_plot_title = "threshold_on_low&medium_vs_high_N_data"
        threshold_plot_labels_by_case={0:"Low&Medium",1:"High"}
        case0_color="cyan"
        linear_regression_threshold2=make_both_ROC_plots(ROC_plot_base_name,
                                 colors2, data, likely_range_for_threshold,
                                 linear_regression_threshold_color2, case0_color,
                                 out_folder, specs, target2,
                                 threshold_plot_title,threshold_plot_labels_by_case)

        #plot the two thresholds together
        threshold_plot_title3 = "threshold_on_low_vs_medium_vs_high_N_data"
        colors3 = ["blue" if m == "High" else "cyan" if m == "Medium" else "gray" for m in true_category]
        threshold_plot_labels_by_case = {0: "Low",1: "Low", 2: "High"}
        fig, ax = plt.subplots(1, 1, figsize=(4, 4))
        alpha = 0.2
        for i in range(0, len(specs)):
            plt.scatter(specs[i], data[i], marker='o', c=colors3[i], alpha=alpha)

        ax.axhline(y=linear_regression_threshold1, color=linear_regression_threshold_color1, linestyle='--',
                   label="thresh1 lr"
                         + " ({0})".format(round(linear_regression_threshold1, 4)))
        ax.axhline(y=linear_regression_threshold2, color=linear_regression_threshold_color2, linestyle='--',
                   label="thresh2 lr"
                         + " ({0})".format(round(linear_regression_threshold2, 4)))
        ax.set(xlabel=" spec time ")
        ax.set(ylabel=" metric ")
        plt.legend(loc=4)
        ax.set(title=threshold_plot_title3)
        plot_file = os.path.join(out_folder, threshold_plot_title3 + ".png")
        plt.savefig(plot_file)
        plt.close()

        return

def make_custom_ROC_plot(X_test, clf, custom_prob_thresholds,
                         file_basename, out_folder, y_test,accuracy_string):

    y_pred_proba=clf.predict_proba(X_test)[:, 1]
    auc = metrics.roc_auc_score(y_test, y_pred_proba)
    fpr_t = [];
    tpr_t = [];
    thresholds_t = [];
    for t in custom_prob_thresholds:
        y_pred_t = (clf.predict_proba(X_test)[:, 1] > t).astype('float')
        result = confusion_matrix(y_test, y_pred_t).ravel()
        tn, fp, fn, tp = result
        total = sum([tn, fp, fn, tp])
        fpr = fp / total
        tpr = tp / total
        fpr_t.append(fpr)
        tpr_t.append(tpr)
        thresholds_t.append(t)
    title = file_basename.replace(".csv", " custom ROC plot")
    fig, ax = plt.subplots(1, 1, figsize=(5, 5))

    plt.plot(fpr_t, tpr_t, label="auc={0},accuracy={1}".format(auc, accuracy_string), marker='x')
    ax.set(xlabel="<-- fpr -->")
    ax.set(ylabel="<-- tpr -->")
    plt.legend(loc=4)
    ax.set(title=title)
    plot_file = os.path.join(out_folder, title.replace(" ", "_") + ".png")
    plt.savefig(plot_file)
    plt.close()

def make_basic_ROC_plot(file_basename, clf, out_folder,
                        X_test, y_test, accuracy_string):

    y_pred_proba = clf.predict_proba(X_test)[::, 1]
    fpr, tpr, thresholds = metrics.roc_curve(y_test, y_pred_proba)
    auc = metrics.roc_auc_score(y_test, y_pred_proba)
    print("thresholds:\t" + str(thresholds))
    title = file_basename.replace(".csv", " ROC plot")
    fig, ax = plt.subplots(1, 1, figsize=(5, 5))
    plt.plot(fpr, tpr, label="auc={0},accuracy={1}".format(auc,accuracy_string), marker='o')
    ax.set(xlabel="<-- fpr -->")
    ax.set(ylabel="<-- tpr -->")
    plt.legend(loc=4)
    ax.set(title=title)
    plot_file = os.path.join(out_folder, title.replace(" ", "_") + ".png")
    plt.savefig(plot_file)
    plt.close()

def make_both_ROC_plots(ROC_plot_base_name, colors, data, likely_range_for_threshold,
                        linear_regression_threshold_color, case0_color, out_folder,
                        specs, target, threshold_plot_title, plot_labels_by_case):
    target_as_array = np.array(target)
    data_as_array = np.array([data]).transpose()
    X_train, X_test, y_train, y_test = train_test_split(data_as_array, target_as_array, test_size=0.33,
                                                        random_state=44)
    clf = LogisticRegression(penalty='l2', C=0.1)
    clf.fit(X_train, y_train)
    y_pred = clf.predict(X_test)
    accuracy = metrics.accuracy_score(y_test, y_pred)
    accuracy_string = str(round(accuracy, 4))
    make_basic_ROC_plot(ROC_plot_base_name, clf, out_folder, X_test, y_test, accuracy_string)
    custom_prob_thresholds = np.arange(0, 2, 0.005)
    make_custom_ROC_plot(X_test, clf, custom_prob_thresholds,
                              ROC_plot_base_name, out_folder, y_test, accuracy_string)
    linear_regression_threshold = get_linear_regession_metric_threshold(clf, likely_range_for_threshold)
    plot_threshold_against_data(colors, data, linear_regression_threshold,
                                     linear_regression_threshold_color, case0_color,
                                     out_folder, threshold_plot_title, specs, plot_labels_by_case)
    print(linear_regression_threshold)
    return linear_regression_threshold

def get_linear_regession_metric_threshold(clf, likely_range_for_threshold):
    model_predictions = clf.predict(np.array([likely_range_for_threshold]).transpose())
    model_prediction_change = [model_predictions[i + 1] - model_predictions[i] for i in
                               range(0, len(likely_range_for_threshold) - 1)]
    i_of_threshold = [i for i in range(0, len(model_prediction_change)) if model_prediction_change[i] == 1]
    linear_regression_threshold = 0.5 * (
                likely_range_for_threshold[i_of_threshold[0]] + likely_range_for_threshold[i_of_threshold[0] + 1])
    return linear_regression_threshold

def plot_threshold_against_data(colors, data, linear_regression_threshold,
                                linear_regression_threshold_color, case0_color,
                                out_folder, plot_title, specs, plot_labels_by_case):
    fig, ax = plt.subplots(1, 1, figsize=(4, 4))
    labeled_0=False
    labeled_1=False
    alpha=0.2
    for i in range(0,len(specs)):

        if colors[i]==case0_color and labeled_0==0:
            plt.scatter(specs[i], data[i], label=plot_labels_by_case[0], marker='o',
                        c=colors[i],alpha=alpha)
            labeled_0=True
        elif colors[i]!=case0_color and labeled_1==0:
            plt.scatter(specs[i], data[i], label=plot_labels_by_case[1], marker='o',
                        c=colors[i],alpha=alpha)
            labeled_1=True
        else:
            plt.scatter(specs[i], data[i], marker='o', c=colors[i],alpha=alpha)


    ax.axhline(y=linear_regression_threshold, color=linear_regression_threshold_color, linestyle='--',
               label="thresh lr"
                     + " ({0})".format(round(linear_regression_threshold, 4)))
    ax.set(xlabel=" spec time ")
    ax.set(ylabel=" metric ")
    plt.legend(loc=4)
    ax.set(title=plot_title)
    plot_file = os.path.join(out_folder, plot_title + ".png")
    plt.savefig(plot_file)
    plt.close()


def get_highN_vs_lowN_truth():

    out_folder = "/home/tamsen/Data/Specks_outout_from_mesx"
    file1="metric3.csv"
    full_path1 = os.path.join(out_folder, file1)

    spc_xs,metric,wgd_sims = read_xs_ys_csv(full_path1)
    low_N_sims= ["Auto","sim37_N0p1", "sim37_N1"]
    med_N_sims=["sim37_N5", "sim35_log"]
    high_N_sims=["sim36_N10","sim37_N20"]

    sims=[]
    specs=[]
    category=[]
    data=[]

    cutoff=60
    for i in range(0,len(wgd_sims)):
        sim = wgd_sims[i]
        spc_time=spc_xs[i]

        if spc_time >= cutoff:
            continue
        if metric[i] <= 0:
            continue
        true_category= categorize_sim(low_N_sims, med_N_sims, high_N_sims, sim)
        sims.append(sim)
        specs.append(spc_time)
        category.append(true_category)
        data.append(math.log(metric[i]))

    return sims,specs,category,data
def get_threshold_based_accuracy(data, target, threshold):
    predictions = [1 if d > threshold else 0 for d in data]
    result = confusion_matrix(target, predictions).ravel()
    tn, fp, fn, tp = result
    print(result)
    accuracy = metrics.accuracy_score(target, predictions)
    print("Accuracy:\t" + str(accuracy))
    return predictions, accuracy, tn, fp, fn, tp

def load_allo_vs_auto_data(csv_file):
    with open(csv_file, "r") as f:
        lines = f.readlines()

    specs = []
    target = []
    data = []
    for l in lines:
        if "truth" in l:
            continue
        if "Allo" in l:
            target.append(1)
        else:
            target.append(0)

        splat = l.split(",")
        specs.append(float(splat[1]))
        data.append(float(splat[2]))

    print("Target:\t" + str(target))
    print("Data:\t" + str(data))
    return specs,data, target


if __name__ == '__main__':
    unittest.main()
