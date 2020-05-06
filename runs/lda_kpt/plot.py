from matplotlib import pyplot as plt
import pandas as pd
import numpy as np
from ast import literal_eval


def converged(data, column='ecut'):
    threshold = 2*(7.349988e-5)
    print("THRESHOLD: {}".format(threshold))
    ref_val = data["energy"].values[-1]

    data["converged"] = np.abs(data["energy"] - ref_val) < threshold
    data["diff"] = np.abs(data["energy"] - ref_val)
    #print(data[["k-points", "converged", "diff"]])
    #print("CONVERGED TO {} meV".format(data["diff"].values[-2]/7.349988e-5))
    return data["diff"].values[-2]/7.349988e-5


def num_atoms(data):
    num_cu = []
    num_sn = []
    num_atom = []
    for i, row in data.iterrows():
        el = row["Element"]
        num_cu.append(1 if((el=="Cu") or (el=="CuSn3")) else 3 \
                      if(el=="Cu3Sn") else 2 if(el=="CuSn") else 0)
        num_sn.append(1 if(el=="Cu3Sn") else 3 \
                      if(el=="CuSn3") else 2 if((el=="Sn") or (el=="CuSn")) else 0)
        num_atom.append(num_cu[-1] + num_sn[-1])

    data["num_cu"] = num_cu
    data["num_sn"] = num_sn
    data["num_atom"] = num_atom
    return data


def formation(data):
    print(len(data))
    data = num_atoms(data)
    cu_ref = data[data["Element"] == "Cu"]["energy"].values[-1]
    sn_ref = data[data["Element"] == "Sn"]["energy"].values[-1]/2

    data["con_formation_energy_per_atom"] = (data["energy"] - data["num_cu"]*cu_ref - \
                                         data["num_sn"]*sn_ref)/data["num_atom"]

    cu_dat = data[data["Element"]=="Cu"][["k-points","energy"]]
    sn_dat = data[data["Element"]=="Sn"][["k-points", "energy"]]
    print(cu_dat)
    print(sn_dat)
    exit(1)

    form_en = []
    for i in range(len(data)):
        kpt = data.iloc[i]["k-points"]
        #print(data.iloc[i]["k-points"])
        if(data.iloc[i]["Element"] == "Cu" or data.iloc[i]["Element"]==0):
            form_en.append(0)
            continue
        try:
            cu_ref_en = cu_dat[cu_dat["k-points"]==kpt]["energy"].values[0]
            sn_ref_en = sn_dat[sn_dat["k-points"]==kpt]["energy"].values[0]
        except IndexError:
            form_en.append(np.nan)
            continue
        #print(data.iloc[i]["k-points"], cu_ref_en, sn_ref_en)
        form_en.append((data.iloc[i]["energy"] - \
                        data.iloc[i]["num_cu"]*cu_ref_en - \
                        data.iloc[i]["num_sn"]*sn_ref_en/2)/data.iloc[i]["num_atom"])
        #print()

    print(form_en)
    data["formation_energy_per_atom"] = form_en
    print(data[["formation_energy_per_atom", "con_formation_energy_per_atom"]])

    #print(data[data["Element"]=="CuSn"][["Element", "k-points", "energy"]])
    #print(data[data["Element"]=="Cu3Sn"][["Element", "k-points", "energy"]])
    #print(data[data["Element"]=="CuSn3"][["Element", "k-points", "energy"]])
    #exit(1)

    #print(data[["Element", "energy", "formation_energy_per_atom"]])
    #print(cu_ref, sn_ref)
    #cusn = data[data["Element"] == 'CuSn']
    #cu3sn = data[data["Element"] == 'Cu3Sn']
    #cusn3 = data[data["Element"] == 'CuSn3']
    #print(cusn[["Element", "energy", "k-points", "formation_energy_per_atom"]])
    #print(cu3sn[["Element", "energy", "k-points", "formation_energy_per_atom"]])
    #print(cusn3[["Element", "energy", "k-points", "formation_energy_per_atom"]])
    #data = num_atoms(data)
    #cu_ref = data[data["Element"] == "Cu"]["energy"].values[-1]
    #sn_ref = data[data["Element"] == "Sn"]["energy"].values[-1]/2

    #data["formation_energy_per_atom"] = (data["energy"] - data["num_cu"]*cu_ref - \
    #                                     data["num_sn"]*sn_ref)/data["num_atom"]

    ##print(data[["Element", "energy", "formation_energy_per_atom"]])
    ##print(cu_ref, sn_ref)
    #cusn = data[data["Element"] == 'CuSn']
    #cu3sn = data[data["Element"] == 'Cu3Sn']
    #cusn3 = data[data["Element"] == 'CuSn3']
    #print(cusn[["Element", "energy", "k-points", "formation_energy_per_atom"]])
    #print(cu3sn[["Element", "energy", "k-points", "formation_energy_per_atom"]])
    #print(cusn3[["Element", "energy", "k-points", "formation_energy_per_atom"]])
    data.to_csv("./final_lda_data.csv")


if __name__ == '__main__':
    data = pd.read_csv("./data_so_far.csv")
    data = pd.concat([data, pd.read_csv("./more_cu_data.csv"),
                      pd.read_csv("./lda_low_kpt.csv")])
    data = data.sort_values('k-points')
    print(data.columns)
    formation(data)
    exit(1)
    #data = pd.concat([data, pd.read_csv("./more_kpt_data.csv"),
    #                  pd.read_csv("./more_cu_kpt_data.csv"),
    #                  pd.read_csv("./alloy_data.csv"),
    #                  pd.read_csv("./more_alloy_data.csv"),
    #                  pd.read_csv("./new_alloy_data.csv"),
    #                  pd.read_csv("./more_new_alloy_data.csv"),
    #                  pd.read_csv("./fill_copper.csv"),
                      #pd.read_csv("./")
    #                  pd.read_csv("./i_hate_copper.csv")])
    #cusn_14 = pd.Series({'Element':'CuSn', 'k-points':14, 'supercell':1,
    #                     'energy':-25299.5787044, 'force':0.0, 'ecut':45})
    #cusn_15 = pd.Series({'Element':'CuSn', 'k-points':15, 'supercell':1,
    #                     'energy':-25299.5823136, 'force':0.0, 'ecut':45})
    
    #cu3sn_14 = pd.Series({'Element':'Cu3Sn', 'k-points':14, 'supercell':1,
    #                     'energy':-25896.3358892, 'force':0.0, 'ecut':45})
    #cu3sn_15 = pd.Series({'Element':'Cu3Sn', 'k-points':15, 'supercell':1,
    #                     'energy':-25896.3330705, 'force':0.0, 'ecut':45})
    #data = data.append(cusn_14, ignore_index=True)
    #data = data.append(cu3sn_14, ignore_index=True)
    #data = data.append(cusn_15, ignore_index=True)
    #data = data.append(cu3sn_15, ignore_index=True)
    #data = data.sort_values('k-points')
    #data.to_csv("./combined_data.csv")
    #print(data)
    #exit(1)

    #data = pd.read_csv("./data.csv")
    #ecut = list(range(40, 71, 5))*2
    #data['ecut'] = pd.Series(ecut)
    #more_data = pd.read_csv("./low_cut_data.csv")
    #data = pd.concat([data, more_data])
    #data = data.sort_values("ecut")
    #print(data.columns)

    for pf in data.groupby("Element"):
        print(pf[0])
        fig, ax = plt.subplots()
        to = converged(pf[1])
        ax.scatter(pf[1]["k-points"], pf[1]["energy"])
        ax.set(title="{0} to {1:.3f} meV".format(pf[0], to))
        plt.savefig("./{}.png".format(pf[0]))
        plt.close()
        #plt.show()
    
