import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
from ast import literal_eval


def most_accurate(data):
    #pass
    bests = []
    for g in data.groupby("Element"):
        #print(g[0], g[1].sort_values("k-points")["con_formation_energy_per_atom"].values[-1])
        #print(g[1].sort_values("k-points").iloc[-1])
        bests.append(g[1].sort_values("k-points").iloc[-1])

    #print(bests[0].index)
    new_dat = pd.DataFrame(bests, columns=bests[0].index)
    print(new_dat[["Element", "con_formation_energy_per_atom"]])
    #best_dat = pd.concat([b for b in bests], axis=1)
    #print(best_dat)
    return new_dat[(new_dat["Element"]=="CuSn") | (new_dat["Element"]=="CuSn3") |
                   (new_dat["Element"]=="Cu3Sn")]


if __name__ == '__main__':
    original = pd.read_csv("./cleaned_bronze_calculations.csv")
    # Get data we can compare to
    #original = original[(original["pretty_formula"] == "Cu") | 
    #              (original["pretty_formula"] == "Cu3Sn") | 
    #              (original["pretty_formula"] == "CuSn") | 
    #              (original["pretty_formula"] == "CuSn3") | 
    #              (original["pretty_formula"] == "Sn")]#[["pretty_formula",
    #                                                   # "formation_energy_per_atom"]]
    sgs = []
    for row in original.iterrows():
        sgs.append(literal_eval(row[1]["spacegroup"])["number"])

    original["sg"] = sgs
    #original = original[
    #            ((original["sg"]==194) & (original["pretty_formula"]=="CuSn")) |
    #            ((original["sg"]==225) & (original["pretty_formula"]=="Cu3Sn")) |
    #            ((original["sg"]==139) & (original["pretty_formula"]=="CuSn3"))
    #]
    print(original[["pretty_formula", "sg", "density", "volume", "num_atoms", "alpha", "beta", "gamma"]])
    exit(1)

    # Attempt to recreate
    recreate_data = pd.read_csv("./runs/recreate/recreated.csv")
    recreate_data["cu_frac"] = recreate_data["num_cu"]/recreate_data["num_atom"]
    print(recreate_data[["Element", "formation_energy_per_atom"]])
    best_recreate_data = most_accurate(recreate_data)
    #exit(1)

    lda_data = pd.read_csv("./runs/lda_kpt/final_lda_data.csv")
    lda_data["cu_frac"] = lda_data["num_cu"]/lda_data["num_atom"]
    best_lda_data = most_accurate(lda_data)

    print(original.columns)
    original["cu_frac"] = original["num_cu"]/original["num_atoms"]

    print(original[["pretty_formula", "formation_energy_per_atom", "not_so_crude"]])
    #exit(1)
    fig, ax = plt.subplots()
    #ax.scatter(best_recreate_data["cu_frac"],
    #           best_recreate_data["con_formation_energy_per_atom"], marker='o', color='b',
    #           edgecolor='k', label="GGA")
    #ax.scatter(best_lda_data["cu_frac"], best_lda_data["con_formation_energy_per_atom"],
    #           marker='h', color='m', edgecolor='k', label="LDA")
    ax.scatter(original["cu_frac"], original["formation_energy_per_atom"], label="Actual",
               marker='s',  color='r', edgecolor='k')
    ax.scatter(original["cu_frac"], original["crude"], label="Crude", marker='p',
               color='g', edgecolor='k')
    ax.scatter(original["cu_frac"], original["not_so_crude"], label="Not so crude", marker='P',
               color='m', edgecolor='k')
    ax.set(xlabel="Cu Fraction", ylabel="Formation Energy per Atom",
           title="Comparing Crude and Not-So-Crude Estimates of Formation Energy")
    ax.legend(loc='best')
    plt.show()

