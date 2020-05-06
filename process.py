import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from ast import literal_eval
import re


def lattice_params(data, kind='bronze'):
    a = []
    b = []
    c = []
    alpha = []
    beta = []
    gamma = []
    for i in data.index:
        d = data.iloc[i]
        sg = literal_eval(d['spacegroup'])['number']

        with open("./{}_cifs/{}_{}.cif".format(kind, d["pretty_formula"], sg), 'r') as fout:
            lines = fout.readlines()
            a.append(float(lines[3].rstrip().split(" ")[-1]))
            b.append(float(lines[4].rstrip().split(" ")[-1]))
            c.append(float(lines[5].rstrip().split(" ")[-1]))
            alpha.append(float(lines[6].rstrip().split(" ")[-1]))
            beta.append(float(lines[7].rstrip().split(" ")[-1]))
            gamma.append(float(lines[8].rstrip().split(" ")[-1]))
        fout.close()

    data["a"] = a
    data["b"] = b
    data["c"] = c
    data["alpha"] = alpha
    data["beta"] = beta
    data["gamma"] = gamma
    return data


def num_atoms(data, kind='bronze'):
    nums = []
    for i in data.index:
        d = data.iloc[i]
        sg = literal_eval(d['spacegroup'])['number']

        with open("./{}_cifs/{}_{}.cif".format(kind, d["pretty_formula"], sg), 'r') as fout:
            last_line = fout.readlines()[-1]
            key = last_line.split(" ")[4]
            nums.append(int(re.split("[A-z]", key)[-1])+1)
        fout.close()
    data["num_atoms"] = nums
    return data

if __name__ == '__main__':
    
    data = pd.read_csv("./bronze_data.csv")
    print(data.columns)
    print(data["band_gap"])
    #exit(1)
    data = lattice_params(data)
    #print(data[["a", "b", "c"]])
    #print(data.columns)
    data = num_atoms(data, 'bronze')
    #print(data[["pretty_formula", "num_atoms", "spacegroup"]])
    #exit(1)
    #data['k-points'] = [2, 3, 5, 3, 1, 3, 5, 5, 5, 3, 2, 5, 3, 1, 2, 5]
    #print(data[["pretty_formula", "num_atoms", "k-points"]])
    #data.to_csv("./bronze_data.csv")
    ncu = []
    nsn = []
    for i in data.index:
        s = data.iloc[i]["pretty_formula"]
        split = re.split("[A-z]", s)

        num_cu = split[2]
        num_cu = 1 if(num_cu == '') else int(num_cu)

        num_sn = split[-1]
        num_sn = 1 if(num_sn == '') else int(num_sn)
        print(num_cu, num_sn)

        base_num = num_sn + num_cu
        ncu.append(int(num_cu/base_num * data.iloc[i]["num_atoms"]))
        nsn.append(int(num_sn/base_num * data.iloc[i]["num_atoms"]))

    print(ncu)
    print(nsn)
    #exit(1)

    data["num_cu"] = ncu
    data["num_sn"] = nsn
    #data.to_csv("./bronze_data.csv")

    new_data = pd.read_csv("./base_data.csv")
    new_en = pd.read_csv("./base_energies.csv")
    #new_data['k-points'] = [5,5,5,5]
    #new_data = num_atoms(new_data, 'base')
    #print(new_data[["pretty_formula", "num_atoms", "k-points"]])
    #new_data.to_csv("./base_data.csv")

    calc = pd.read_csv("./calculated_energies.csv")
    #print(data.columns)
    #print(calc.columns)
    #print(calc[["pretty_formula", "energy"]])

    #joined = pd.concat([data, calc], axis='columns')
    #print(joined.columns)
    #print(joined[["pretty_formula", "energy", "formation_energy_per_atom"]])

    #print(new_data.columns)
    cu_ref = new_en[new_en["pretty_formula"]=='Cu']["energy"].values
    sn_ref = new_en[new_en["pretty_formula"]=='Sn']["energy"].values/2

    ncu_ref = -1592.724318
    nsn_ref = -4048.554436/2

    print(cu_ref, sn_ref)
    #exit(1)
    #print(cu_ref, sn_ref)

    #print(data[["pretty_formula", "num_atoms", "num_cu", "num_sn"]])
    #print(calc[["pretty_formula", "energy"]])
    #print(calc)
    print(calc[["pretty_formula", "energy"]])
    print(data[["pretty_formula", "num_cu", "num_sn", "num_atoms"]])
    print(cu_ref, sn_ref)
    calc["crude"] = (calc["energy"] - data["num_cu"]*cu_ref - data["num_sn"]*sn_ref)/data["num_atoms"]
    calc["not_so_crude"] = (calc["energy"] - data["num_cu"]*ncu_ref - data["num_sn"]*nsn_ref)/data["num_atoms"]
    data = pd.concat([data, calc[calc.columns]], axis=1)
    print(data.columns)
    #exit(1)
    print(data[["pretty_formula", "formation_energy_per_atom", "crude", "not_so_crude"]])
    #print(calc)
    data.to_csv("./cleaned_bronze_calculations.csv")

