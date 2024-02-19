from matplotlib import pyplot as plt
import matplotlib.colors as mcolors


def get_Ks_from_file(cvs_file):
    ks_results = []
    with open(cvs_file, "r") as f:
        lines = f.readlines()

    for l in lines[1:len(lines)]:
        data_splat = l.split(",")
        # print(data_splat)
        ks = float(data_splat[2])
        ks_results.append(ks)

    return ks_results


def histogram_ks(ks_data, out_file):
    ks_reduced = [ks / 1.0 for ks in ks_data if ks < 10]
    fig = plt.figure(figsize=(10, 10), dpi=100)
    n_bins = 50
    n, bins, patches = plt.hist(ks_reduced, bins=n_bins, facecolor='b', alpha=0.25, label='histogram data')
    plt.axvline(x=2, color='r', linestyle='--', label="WGD")
    plt.axvline(x=3, color='b', linestyle='--', label="SPEC")
    plt.legend()
    plt.xlim(0, 11)
    plt.xlabel("Node dist (MYA)")
    plt.ylabel("Num nodes")
    plt.title("Node-distance histogram for " + "polyploid.species_name")
    plt.savefig(out_file)
    plt.clf()
    plt.close()

def csv_to_histogram(f1):
    ks=get_Ks_from_file(f1)
    hist_file=f1.replace(".csv",".png")
    histogram_ks(ks,hist_file)