import numpy, os
import numpy as np
import matplotlib.pyplot as plt
from labutil.plugins.pwscf import run_qe_pwscf, PWscf_inparam, parse_qe_pwscf_output
from labutil.objects import Struc, Dir, ase2struc, Kpoints, PseudoPotential
from ase.spacegroup import crystal
from ase.io import write, read
from ase.build import bulk, make_supercell
import pandas as pd
from ast import literal_eval
import re
import time


#def make_struc(alat, cell='hcp', mag='non'):
def make_struc(formula, spacegroup, sc=1):
    """
    Creates the crystal structure using ASE.
    :param alat: Lattice parameter in angstrom
    :return: structure object converted from ase
    """
    #src_dir = 'bronze' if('Cu' in formula) else 'steel'
    src_dir = 'bronze' if(('C' in formula) and ('S' in formula)) else 'base'
    cif_in = read("./{}_cifs/{}_{}.cif".format(src_dir, formula, spacegroup))
    supercell = make_supercell(cif_in, [[sc,0,0],[0,sc,0],[0,0,sc]])
    structure = Struc(ase2struc(supercell))
    return structure


def compute_energy(nk, ecut, elements, formula, spacegroup, sc=1,
                   directory="final_project", ncpu=8, functional="GGA"):
    """
    Make an input template and select potential and structure, and the path where to run
    """

    # Base Metal
    if('Cu' in elements):
        e1 = 'Cu'
        if(functional=="GGA"):
            pot1 = 'Cu.rel-pbe-spn-kjpaw_psl.1.0.0.UPF'
        elif(functional=="LDA"):
            pot1 = 'Cu.pz-dn-rrkjus_psl.0.2.UPF'
            #pot1 = 'Cu.rel-pbe-spn-kjpaw_psl.1.0.0.UPF'

    # Alloying agent
    if('Sn' in elements):
        e2 = 'Sn'
        if(functional=="GGA"):
            pot2 = 'Sn.rel-pbe-dn-kjpaw_psl.1.0.0.UPF'
        elif(functional=="LDA"):
            pot2 = 'Sn.pz-dn-rrkjus_psl.0.2.UPF'
        if('Cu' not in elements):
            e1 = 'Sn'
            if(functional=="GGA"):
                pot1 = 'Sn.rel-pbe-dn-kjpaw_psl.1.0.0.UPF'
            elif(functional=="LDA"):
                pot1 = 'Sn.pz-dn-rrkjus_psl.0.2.UPF'


    # Put it together
    pot1path = './pseudopot/' + pot1
    if(len(elements) == 1):
        pseudopots = {e1: PseudoPotential(path=pot1path, ptype='paw', element=e1,
                                        functional='GGA', name=pot1)
                  }
    elif(len(elements) == 2):
        pot2path = './pseudopot/' + pot2
        pseudopots = {e1: PseudoPotential(path=pot1path, ptype='paw', element=e1,
                                        functional='GGA', name=pot1),
                  e2: PseudoPotential(path=pot2path, ptype='paw', element=e2,
                                        functional='GGA', name=pot2)
                  }

    # Get structure
    struc = make_struc(formula, spacegroup, sc)
    kpts = Kpoints(gridsize=[nk, nk, nk], option='automatic', offset=False)

    # Directory
    dirname = '{}_{}_{}_{}_{}'.format(formula, spacegroup, nk, sc, ecut)

    # Path
    runpath = Dir(path=os.path.join(os.environ['WORKDIR'], "final_project/"+directory,
                  dirname))

    # Create input
    input_params = PWscf_inparam({
        'CONTROL': {
            'calculation': 'scf',
            #'pseudo_dir': os.environ['QE_POTENTIALS'],
            'pseudo_dir': '/home/bond/Work/final_project/pseudopot',
            'outdir': runpath.path,
            'tstress': True,
            'tprnfor': True,
            'disk_io': 'none',
        },
        'SYSTEM': {
            'ecutwfc': ecut,
            'ecutrho': ecut * 10,
            'nspin': 4,
            'starting_magnetization(1)': 0,
            'starting_magnetization(2)': 0,
            'occupations': 'smearing',
            'smearing': 'mp',
            'degauss': 0.02,
            'ntyp': 2,
            'lspinorb': True
             },
        'ELECTRONS': {
            'diagonalization': 'david',
            'mixing_beta': 0.5,
            'conv_thr': 1e-7,
        },
        'IONS': {},
        'CELL': {},
        })

    # Run and parse
    #print(struc)
    #exit(1)
    output_file = run_qe_pwscf(runpath=runpath, struc=struc,  pseudopots=pseudopots,
                               params=input_params, kpoints=kpts, ncpu=ncpu)
    output = parse_qe_pwscf_output(outfile=output_file, col=True)
    return output


def num_atoms(data):

    nums = []
    for i in bronze_data.index:
        bd = bronze_data.iloc[i]
        sg = literal_eval(bd['spacegroup'])['number']
        cif_in = read("./{}_cifs/{}_{}.cif".format("bronze",
                                                   bd["pretty_formula"], sg))
        with open("./{}_cifs/{}_{}.cif".format("bronze",
                         bd["pretty_formula"], sg), 'r') as fout:
            last_line = fout.readlines()[-1]
            key = last_line.split(" ")[4]
            nums.append(int(re.split("[A-z]", key)[-1]))
        fout.close()

    #print(nums)
    data["num_atoms"] = nums
    #print(data[["pretty_formula", "num_atoms"]])

    return data


def convergence():
    #ks = [5,6,7,8,9,10,11,12,13]
    #ks = [14,15,16]
    #ks = [2,3,4]
    ks = [5]
    sc = [1]#,2,3]
    #element = ['Cu', 'Sn']#, 'CuSn']
    #element = ['Cu']
    element = ['Cu', 'Sn', 'CuSn', 'Cu3Sn', 'CuSn3']
    element = ['Sn', 'Cu']
    ecuts = [45]#,45,50,55,60,65,70]
    funct = 'GGA'

    kpt = []
    supercell = []
    compound = []
    energy = []
    pressure = []
    force = []
    cutoff = []

    for s in sc:
        for k in ks:
            #if((s==1 and k<6) or (s==2 and (k<8 and k>2)) or (s==3 and k>3)):
            #    continue
            for e in element:
                for ecut in ecuts:
                    print("{}: {} kpts, {} supercell {} ecut".format(e, k, s, ecut))
                    sg = 225 if(e=='Cu') else 227 if(e=='Sn') else 194 if(e=='CuSn') \
                            else 139 if(e=='CuSn3') else 225 #if(e=='Cu3Sn')
                    #formula = 'CuSn' if(len(e)==2) else e
                    el = [e] if(len(e)==2) else ['Cu', 'Sn']
                    start = time.time()
                    output = compute_energy(nk=k, ecut=ecut, elements=el, formula=e,
                                    spacegroup=sg, directory='lda_kpt', sc=s, ncpu=4,
                                    functional=funct)

                    print("TIME: {}".format(time.time() - start))
                    kpt.append(k)
                    supercell.append(s)
                    compound.append(e)
                    energy.append(output["energy"])
                    #pressure.append(output["pressure"])
                    force.append(output["force"])
                    cutoff.append(ecut)
                d = pd.DataFrame({
                        'Element': compound,
                        'k-points': kpt,
                        'supercell': supercell,
                        'energy': energy,
                        #'pressure': pressure,
                        'force': force,
                        'ecut': cutoff
                })
    
                #d.to_csv("./lda_kpt/small.csv")
                d.to_csv("./recreate/small.csv")
            print()
            #break
        #break

    d = pd.DataFrame({
            'Element': compound,
            'k-points': kpt,
            'supercell': supercell,
            'energy': energy,
            #'pressure': pressure,
            'force': force,
            'ecut': cutoff
    })

    d.to_csv("./recreate/just_sn_5_kpt.csv")


if __name__ == '__main__':
    bronze_data = pd.read_csv("./bronze_data.csv")
    bronze_data = num_atoms(bronze_data)
    sg_nums = []
    for i in bronze_data.index:
        sg_nums.append(literal_eval(bronze_data.loc[i]['spacegroup'])['number'])
    bronze_data["sg_num"] = sg_nums
    #print(bronze_data[["pretty_formula", "num_atoms", "sg_num"]])
    #exit(1)
    convergence()
    exit(1)

    #cif_in = read("./steel_cifs/Fe4C_87.cif")
    #struc = ase2struc(cif_in)

    bronze_data = pd.read_csv("./base_data.csv")
    bronze_data = num_atoms(bronze_data)

    #bd = bronze_data.iloc[2]
    #sg = literal_eval(bd['spacegroup'])
    #print(bd)
    #print(sg['number'])
    #exit(1)
    calc_en = []
    pf = []
    force = []
    pressure = []
    for i in bronze_data.index:

        bd = bronze_data.iloc[i]
        #print(bd["pretty_formula"])
        sg = literal_eval(bd['spacegroup'])
        nk = 5 if(bd['num_atoms'] < 10) else 3 if(bd['num_atoms'] < 30) else 2 \
               if(bd['num_atoms'] < 40) else 1
        output = compute_energy(nk=2, ecut=40, elements=bd['elements'],
                            formula=bd['pretty_formula'], spacegroup=sg['number'],
                            directory='base_calculations')

        pf.append(bd["pretty_formula"])
        calc_en.append(output["energy"])
        force.append(output["force"])
        pressure.append(output["pressure"])

    nd = pd.DataFrame({'pretty_formula': pf, 'energy': calc_en,
                       'force':force, 'pressure':pressure})
    print(nd)
    nd.to_csv("./base_energies.csv")
    
    #print(bronze_data["magnetization"])
    # put here the function that you actually want to run
    #lattice_scan()
    #k_point_scan()
    #magnetic_gs()
