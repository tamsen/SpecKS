import math
import os
import unittest

import numpy as np
from sklearn import svm, datasets
#import scikitplot as skplt
import matplotlib.pyplot as plt
from sklearn import metrics
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import confusion_matrix
from sklearn.model_selection import train_test_split
from sklearn.datasets import load_breast_cancer
import matplotlib.pyplot as plt

from results_viewer.allo_auto_predictor import read_xs_ys_csv, categorize_sim


#https://stackoverflow.com/questions/25009284/how-to-plot-roc-curve-in-python
# https://stackoverflow.com/questions/28716241/controlling-the-threshold-in-logistic-regression-in-scikit-learn

class TestLinearRegression(unittest.TestCase):
    def test_linear_regression_on_allo_vs_auto_data(self):

        out_folder = "/home/tamsen/Data/Specks_outout_from_mesx"
        file_basename= "allopolyploid_index_truth_vs_predictions.csv"
        csv_file=os.path.join(out_folder,file_basename)
        data, target = load_allo_vs_auto_data(csv_file)

        #threshold=0.39
        threshold = 0.9284
        predictions, accuracy, tn, fp, fn, tp = get_threshold_based_accuracy(data, target, threshold)

        target_as_array=np.array(target)
        data_as_array=np.array([data]).transpose()
        X_train, X_test, y_train, y_test = train_test_split(data_as_array, target_as_array, test_size=0.33, random_state=44)
        clf = LogisticRegression(penalty='l2', C=0.1)
        clf.fit(X_train, y_train)
        y_pred = clf.predict(X_test)
        accuracy = metrics.accuracy_score(y_test, y_pred)
        print("Lin Regression Accuracy", str( accuracy))
        print("Lin Regression Coeff",str(clf.coef_))
        #Lin Regression Coeff [[0.92845063]]
        #
        y_pred_proba = clf.predict_proba(X_test)[::, 1]
        fpr, tpr, thresholds = metrics.roc_curve(y_test, y_pred_proba)
        auc = metrics.roc_auc_score(y_test, y_pred_proba)
        print("thresholds:\t" + str(thresholds))
        title=file_basename.replace(".csv"," ROC plot")
        fig, ax = plt.subplots(1, 1, figsize=(5, 5))
        plt.plot(fpr, tpr, label="data 1, auc=" + str(auc),  marker='o')
        ax.set(xlabel="<-- fpr -->")
        ax.set(ylabel="<-- tpr -->")
        plt.legend(loc=4)
        ax.set(title=title)
        plot_file = os.path.join(out_folder, title.replace(" ","_")+".png")
        plt.savefig(plot_file)
        plt.close()

        fpr_t=[]; tpr_t=[]; thresholds_t=[];
        custom_thresholds = np.arange(0, 2, 0.005)
        for t in custom_thresholds:
            y_pred_t = (clf.predict_proba(X_test)[:, 1] > t).astype('float')
            result = confusion_matrix(y_test, y_pred_t).ravel()
            tn, fp, fn, tp = result
            total=sum([tn,fp,fn,tp])
            fpr=fp/total
            tpr=tp/total
            fpr_t.append(fpr)
            tpr_t.append(tpr)
            thresholds_t.append(t)

        title = file_basename.replace(".csv", " custom ROC plot")
        fig, ax = plt.subplots(1, 1, figsize=(5, 5))
        plt.plot(fpr_t, tpr_t, label="data 1 by t",marker="x")
        ax.set(xlabel="<-- fpr -->")
        ax.set(ylabel="<-- tpr -->")
        plt.legend(loc=4)
        ax.set(title=title)
        plot_file = os.path.join(out_folder, title.replace(" ", "_") + ".png")
        plt.savefig(plot_file)
        plt.close()

    def test_linear_regression_on_high_vs_low_N_data(self):

        out_folder = "/home/tamsen/Data/Specks_outout_from_mesx"
        sims,specs,true_category,data=get_highN_vs_lowN_truth()

        #for i in range(0,len(sims)):
        #    print([sims[i],specs[i],true_category[i],data[i]])

        target = [0 if m == "Low"
                  else 1 for m in true_category]

        threshold=0.01
        predictions = [1 if d > threshold else 0 for d in data]
        result = confusion_matrix(target, predictions).ravel()
        tn, fp, fn, tp = result
        print(result)
        accuracy = metrics.accuracy_score(target, predictions)
        print("Accuracy:\t" + str(accuracy))

        #plot it
        #colors = ["gray" if m == "Low"
        #          else "cyan" if m == "Medium"
        #            else "blue" for m in true_category]

        colors = ["gray" if m == 0
                  else "cyan" for m in target]

        #LinearRegression fits a linear model with coefficients w = (w1, …, wp) to
        # minimize the residual sum of squares between the observed targets in the dataset,
        # and the targets predicted by the linear approximation.

        target_as_array = np.array(target)
        data_as_array = np.array([data]).transpose()
        X_train, X_test, y_train, y_test = train_test_split(data_as_array, target_as_array, test_size=0.33,
                                                            random_state=44)
        clf = LogisticRegression(penalty='l2', C=0.1)
        clf.fit(X_train, y_train)
        y_pred = clf.predict(X_test)
        accuracy = metrics.accuracy_score(y_test, y_pred)
        print("Lin Regression Accuracy", str(accuracy))
        print("Lin Regression Coeff", str(clf.coef_))
        print("Lin Regression Int", str(clf.intercept_))
        #return
        test_m = np.arange(-10, 0, 1)
        probs=clf.predict(np.array([test_m]).transpose())
        for i in range(0,len(test_m)):
            #print("metric {0} has p {1} prob of being HIGH or MEDIUM".format(test_m[i],probs[i][1]))
            print("metric {0} has p {1} prob of being HIGH or MEDIUM".format(test_m[i],probs[i]))
        print("clf params:" + str(clf.get_params()))

        fig, ax = plt.subplots(1, 1, figsize=(4, 4))
        plt.scatter(specs, data, label='label', marker='o', c=colors)

        ax.axhline(y=threshold, color='cyan', linestyle='--', label="lvm disc. criteria"
                                                                                          + " ({0})".format(
            round(threshold, 2)))
        ax.set(xlabel=" spec time ")
        ax.set(ylabel=" metric ")
        plt.legend(loc=4)
        ax.set(title='title')
        plot_file = os.path.join(out_folder, "metric3_catagories2.png")
        plt.savefig(plot_file)
        plt.close()

        return predictions, accuracy, tn, fp, fn, tp



        #catagory_0 = ["Low","Medium"]
        #catagory_1 =["High"]
        #optimal_threshold_Medium_vs_High = perform_linear_regression(
        #    catagory_0, catagory_1, true_category,data, out_folder)


        #colors = ["gray" if m < optimal_threshold_Low_vs_Medium
        #        else "blue" if  m > optimal_threshold_Medium_vs_High
        #        else "cyan" for m in data]

        colors = ["gray" if m == "Low"
                  else "cyan" if m == "Medium"
                    else "blue" for m in true_category]



        fig, ax = plt.subplots(1, 1, figsize=(4, 4))
        plt.scatter(specs,data, label='label', marker='o', c=colors)

        ax.axhline(y=optimal_threshold_Low_vs_Medium, color='cyan', linestyle='--', label="lvm disc. criteria"
                                                                                       + " ({0})".format(
            round(optimal_threshold_Low_vs_Medium, 2)))

        #ax.axhline(y=optimal_threshold_Medium_vs_High , color='purple', linestyle='--', label="mvh disc. criteria"
        #                                                                                  + " ({0})".format(
        #    round(optimal_threshold_Medium_vs_High , 2)))

        ax.set(xlabel=" spec time ")
        ax.set(ylabel=" metric ")
        plt.legend(loc=4)
        ax.set(title='title')
        plot_file = os.path.join(out_folder, "metric3_with_thresholds.png")
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
def perform_linear_regression(catagory_0, catagory_1, true_category, raw_data, out_folder):

    test_name="".join(catagory_0) +"vs" + "".join(catagory_1)
    data, target = load_high_vs_low_data(true_category, raw_data, catagory_0, catagory_1)
    target_as_array = np.array(target)
    data_as_array = np.array([data]).transpose()
    X_train, X_test, y_train, y_test = train_test_split(data_as_array, target_as_array, test_size=0.33,
                                                        random_state=44)
    clf = LogisticRegression(penalty='l2', C=0.1)
    clf.fit(X_train, y_train)
    y_pred = clf.predict(X_test)
    accuracy = metrics.accuracy_score(y_test, y_pred)
    print("Lin Regression Accuracy", str(accuracy))
    print("Lin Regression Coeff", str(clf.coef_))

    y_pred_proba = clf.predict_proba(X_test)[::, 1]
    fpr, tpr, thresholds = metrics.roc_curve(y_test, y_pred_proba)
    auc = metrics.roc_auc_score(y_test, y_pred_proba)
    print("thresholds:\t" + str(thresholds))
    title = test_name + "N ROC plot"
    fig, ax = plt.subplots(1, 1, figsize=(4, 4))
    label = "data 1, auc=" + str(round(auc,3)) + "\nopt. threshold: " + str(round(clf.coef_[0][0],3))
    plt.plot(fpr, tpr, label=label, marker='o')
    ax.set(xlabel="<-- fpr -->")
    ax.set(ylabel="<-- tpr -->")
    plt.legend(loc=4)
    ax.set(title=title)
    plot_file = os.path.join(out_folder, title.replace(" ", "_") + "_"+test_name + ".png")
    plt.savefig(plot_file)
    plt.close()

    optimal_threshold=clf.coef_[0][0]
    return optimal_threshold

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
        data.append(float(splat[2]))

    print("Target:\t" + str(target))
    print("Data:\t" + str(data))
    return data, target

def load_high_vs_low_data(true_category, raw_data,catagory_0,catagory_1):

    target = []
    data = []
    for i in range(0,len(true_category)):

        if true_category[i] in catagory_0:
            target.append(0)
        else:
            target.append(1)

        data.append(float(raw_data[i]))

    print("Target:\t" + str(target))
    print("Data:\t" + str(data))
    return data, target


if __name__ == '__main__':
    unittest.main()
