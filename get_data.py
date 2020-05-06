from pymatgen.ext.matproj import MPRester
from matplotlib import pyplot as plt
import ast
import pandas as pd
import numpy as np


def get_data(key):
    with MPRester(key) as m:
        compounds = m.query(criteria={'elements':{#"$in":['Cu', 'Sn', 'Fe', 'C', 'Al', 'Ag'],
                                                  '$all':['Cu']},
                                      'nelements':{"$eq":2},
                                      'e_above_hull':{"$eq":0.0}},
                            properties=['elements', 'pretty_formula', 'e_above_hull',
                                        'formation_energy_per_atom', 'spacegroup', 'density',
                                        'nsites', 'cif', 'volume', 'band_gap',
                                        'total_magnetization', 'elasticity', 'piezo', 'diel'])

    data = pd.DataFrame(compounds, columns=['pretty_formula', 'elements', 'e_above_hull',
                                            'formation_energy_per_atom', 'spacegroup',
                                            'density', 'nsites', 'cif', 'volume', 'band_gap',
                                            'total_magnetization', 'elasticity', 'piezo',
                                            'diel'])
    data.to_csv("./cu_all_data.csv", index=False)

if __name__ == '__main__':
    #key = str(open("./API", 'r').read().split()[1])
    #get_data(key)

    data = pd.read_csv("./cu_all_data.csv")
    #print(data[["pretty_formula", "formation_energy_per_atom"]])

    #print(np.unique(data["elements"]))
    #print(data[["pretty_formula", "nsites", "formation_energy_per_atom", "cif"]])
    #print(data.loc[0]["cif"])
    
    for i in data.index:
        #space_group = ast.literal_eval(data.iloc[i]["spacegroup"])
        #file_name = data.iloc[i]["pretty_formula"]
        #with open("./cu_all_compounds/{}_{}.cif".format(file_name, space_group['number']), 'w') \
        #  as fout:
            #fout.write(data.iloc[i]["cif"])
        #fout.close()
        print(data.iloc[i]["cif"].split("\n")[-2].split(" ")[4])

    fig, ax = plt.subplots()

