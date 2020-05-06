from matplotlib import pyplot as plt
import pandas as pd
import numpy as np


def converged(data, column='ecut'):
    threshold = 2*(7.349988e-5)
    ref_val = data["energy"].values[-1]
    data["converged"] = np.abs(data["energy"] - ref_val) < threshold
    data["diff"] = np.abs(data["energy"] - ref_val)
    print(data[["k-points", "converged", "diff"]])


if __name__ == '__main__':
    #data = pd.read_csv("./kpt_data.csv")
    #data = pd.concat([data, pd.read_csv("./more_kpt_data.csv")])
    #data = data.sort_values('k-points')

    data = pd.read_csv("./data.csv")
    ecut = list(range(40, 71, 5))*2
    data['ecut'] = pd.Series(ecut)
    more_data = pd.read_csv("./low_cut_data.csv")
    data = pd.concat([data, more_data])
    data = data.sort_values("ecut")
    print(data.columns)

    fig, ax = plt.subplots()
    e_data = data[data["Element"]=="Cu"]
    converged(e_data)
    ax.scatter(e_data["ecut"], e_data["energy"])
    plt.show()
    pass
