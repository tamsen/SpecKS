import os
import unittest
from sklearn import svm, datasets
#import scikitplot as skplt
import matplotlib.pyplot as plt
from sklearn import metrics
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import confusion_matrix
from sklearn.model_selection import train_test_split
from sklearn.datasets import load_breast_cancer
import matplotlib.pyplot as plt
#https://stackoverflow.com/questions/25009284/how-to-plot-roc-curve-in-python
# https://stackoverflow.com/questions/28716241/controlling-the-threshold-in-logistic-regression-in-scikit-learn

class TestLinearRegression(unittest.TestCase):
    def test_linear_regression(self):

        out_folder = "/home/tamsen/Data/Specks_outout_from_mesx"
        breast_cancer = load_breast_cancer()
        X = breast_cancer.data
        y = breast_cancer.target
        X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.33, random_state=44)
        clf = LogisticRegression(penalty='l2', C=0.1)
        clf.fit(X_train, y_train)
        y_pred = clf.predict(X_test)

        #Accuracy= the fraction of correctly classified samples (float)
        print("Accuracy", metrics.accuracy_score(y_test, y_pred))
        #Accuracy = TP + TN / TOTAL

        y_pred_proba = clf.predict_proba(X_test)[::, 1]
        fpr, tpr, thresholds = metrics.roc_curve(y_test, y_pred_proba)
        auc = metrics.roc_auc_score(y_test, y_pred_proba)
        #print("thresholds:" + str(thresholds))
        print("\n\n")
        print("\tthres\tfpr\ttpr")
        for i in range(0,len(thresholds)):
            #accuracy =
            thresh_str=str(round(thresholds[i],2))
            fpr_str=str(round(thresholds[i],2))
            tpr_str=str(round(thresholds[i],2))
            print("\t{0}\t{1}\t{2}".format(thresh_str,fpr_str,tpr_str))

        #the confusion matrix for a given threshold
        threshold = thresholds[0]
        y_pred_t = (clf.predict_proba(X_test)[:, 1] > threshold).astype('float')
        result=confusion_matrix(y_test, y_pred_t).ravel()
        tn, fp, fn, tp =result
        print(result)

        fig, ax = plt.subplots(1, 1, figsize=(5, 5))
        plt.plot(fpr, tpr, label="data 1, auc=" + str(auc))
        ax.set(xlabel="<-- fpr -->")
        ax.set(ylabel="<-- tpr -->")
        plt.legend(loc=4)
        #plt.show()
        plot_file = os.path.join(out_folder, "ExampleROCplot.png")
        plt.savefig(plot_file)
        plt.close()


    def test_linear_regression2(self):
        breast_cancer = load_breast_cancer()

        X = breast_cancer.data
        y = breast_cancer.target
        X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.33, random_state=44)
        clf = LogisticRegression(penalty='l2', C=0.1)
        clf.fit(X_train, y_train)
        y_pred = clf.predict(X_test)
        print("Accuracy", metrics.accuracy_score(y_test, y_pred))
        y_probas= clf.predict_proba(X_test)[::, 1]

        #skplt.metrics.plot_roc_curve(y, y_probas)
        #plt.show()


if __name__ == '__main__':
    unittest.main()
