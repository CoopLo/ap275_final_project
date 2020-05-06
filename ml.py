import pandas as pd
#import numpy as np
from matplotlib import pyplot as plt
from sklearn.neural_network import MLPRegressor
from sklearn.kernel_ridge import KernelRidge as KR
from sklearn.linear_model import Lasso
from sklearn.model_selection import LeaveOneOut, cross_val_score
from ast import literal_eval
from itertools import combinations
#from numba import jit
import time
import autograd as ad
from autograd import numpy as np


def leave_one_out(data, model, cols):
    target = 'formation_energy_per_atom'
    loo = LeaveOneOut()
    rmse = []
    cv_scores = []
    start = time.time()
    for train_index, test_index in loo.split(data[cols]):

        # Split dataset
        x_train, x_test = data[cols].loc[train_index], data[cols].loc[test_index]
        y_train, y_test = data[target].loc[train_index], data[target].loc[test_index]

        # Fit model
        model.fit(x_train, y_train)

        # Get metrics
        pred_y = model.predict(x_test)
        rmse.append((pred_y[0] - y_test.values[0])**2)
        #cv_scores.append(cross_val_score(model, x_test, y_test, cv=1))

    print("\nAVERAGE RMSE: {}, IN {} SECONDS".format(np.sqrt(np.average(rmse)),
                                                     time.time()-start))
    #print(np.average(cv_scores))
    #exit(1)
    return np.sqrt(np.average(rmse))
    #print(y_test.values.reshape(-1,1), pred_y.reshape(-1,1))
    #return model.score(y_test.values.reshape(-1,1), pred_y.reshape(-1,1))


def nn(data, cols=["num_cu", "num_sn", "crude"]):
    model = MLPRegressor(activation='tanh', solver='adam', learning_rate='adaptive',
                 hidden_layer_sizes=(100,100), verbose=False,
                 max_iter=10000, tol=1e-12)

    #model.fit(data[cols], data["formation_energy_per_atom"])
    #score = model.score(data[cols], data["formation_energy_per_atom"])
    #print("NN SCORE: {}".format(score))
    #print("NN SCORE: {}".format(model.score(data[cols], data["formation_energy_per_atom"])))

    #pred_y = model.predict(data[cols])
    #fig, ax = plt.subplots()
    #ax.scatter(data["num_cu"], data["formation_energy_per_atom"], label="True")
    #ax.scatter(data["num_cu"], pred_y, label="Predicted")
    #ax.legend(loc='best')
    #plt.show()
    return leave_one_out(data, model, cols)
    #return score


def kernel_regression(data, cols=["crude"]):
    sigma = 0.1
    def kernel(x1, x2):
        return np.exp(-np.linalg.norm(x1-x2)**2/(2*sigma**2))

    model = KR(kernel='rbf')
    #model.fit(data[cols], data["formation_energy_per_atom"])
    #score = model.score(data[cols], data["formation_energy_per_atom"])
    #print("KRR SCORE: {}".format(score))
    #print("KRR SCORE: {}".format(model.score(data[cols], data["formation_energy_per_atom"])))

    #pred_y = model.predict(data[cols])
    #fig, ax = plt.subplots()
    #ax.scatter(data["num_cu"], data["formation_energy_per_atom"], label="True", alpha=0.5)
    #ax.scatter(data["num_cu"], pred_y, label="Predicted", alpha=0.5)
    #ax.legend(loc='best')
    #plt.show()
    return leave_one_out(data, model, cols)
    #return score


def lasso(data, cols=["crude"]):
    model = Lasso()
    #model.fit(data[cols], data["formation_energy_per_atom"])
    #score = model.score(data[cols], data["formation_energy_per_atom"])
    #print("LASSO SCORE: {}".format(score))
    #print("LASSO SCORE: {}".format(model.score(data[cols], data["formation_energy_per_atom"])))

    #pred_y = model.predict(data[cols])
    #fig, ax = plt.subplots()
    #ax.scatter(data["crude"], data["formation_energy_per_atom"], label="True", alpha=0.5)
    #ax.scatter(data["crude"], pred_y, label="Predicted", alpha=0.5)
    #ax.legend(loc='best')
    #plt.show()
    return leave_one_out(data, model, cols)
    #return score


if __name__ == '__main__':
    data = pd.read_csv("./cleaned_bronze_calculations.csv")
    #print(data[["pretty_formula", "spacegroup", "total_magnetization"]])
    #exit(1)
    #print(data.columns)
    #exit(1)
    data["diff"] = data["formation_energy_per_atom"] - data["crude"]

    #fig, ax = plt.subplots()
    #ax.scatter(data["num_cu"], data["formation_energy_per_atom"],
    #           label="High Quality Calculations", edgecolor='k', marker='p', color='#e15928')

    #ax.scatter(data["num_cu"], data["crude"], label="Crude Estimate", edgecolor='k',
    #           marker='s', color='#7570b3')

    #ax.set(xlabel="Number of Cu Atoms in Primitive Cell",
    #       ylabel="Formation Energy per Atom (eV)",
    #       title="Crude Estimate vs. High Quality Calculations of Cu-Sn Compounds")
    #ax.legend(loc='best')
    #plt.show()
    #print(literal_eval(data.iloc[0]["spacegroup"]))
    #print(data.iloc[0]['elasticity'])

    cols = ['density', 'volume', 'total_magnetization',
            'num_cu', 'num_sn', 'lattice', 'angles']
    target = 'formation_energy_per_atom'

    best_score = 99999
    best_col_combination = None
    best_model = None

    #loo = LeaveOneOut()
    #rmse = []
    #for train_index, test_index in loo.split(data[cols]):
    #    print("TRAIN: {}, TEST: {}".format(train_index, test_index))
    #    x_train, x_test = data[cols].loc[train_index], data[cols].loc[test_index]
    #    y_train, y_test = data[target].loc[train_index], data[target].loc[test_index]

    #    model.fit(x_train, y_train)
    #    pred_y = model.predict(x_test)
    #    rmse.append((pred_y[0] - y_test.values[0])**2)
    #print("\nAVERAGE RMSE: {}".format(np.average(rmse)))
    #exit(1)

    nn_scores = []
    krr_scores = []
    lasso_scores = []
    nc_nn_scores = []
    nc_krr_scores = []
    nc_lasso_scores = []

    col_combs = []
    for i in range(1, len(cols)+1):
        if(i > 1):
            continue
        for c in list(combinations(cols, i)):
            col_combs.append(c)

            # Sub in lattice parameters
            if('lattice' in c):
                c = np.array(c)
                c = c[c!='lattice']
                c = list(c)
                c.append('a')
                c.append('b')
                c.append('c')

            # Sub in lattice angles
            if("angles" in c):
                c = np.array(c)
                c = c[c!='angles']
                c = list(c)
                c.append('alpha')
                c.append('beta')
                c.append('gamma')

            c = list(c)

            no_crude_c = np.copy(c)
            c.append("crude")

            #print("COLUMNS: {}".format(c))
            print("COLUMNS: {}".format(no_crude_c))
            print("WITH APPROXIMATION")
            avg_nnc = 0
            for i in range(10):
                avg_nnc += 100
                #avg_nnc += nn(data, cols=c)/10
            krrc = kernel_regression(data, cols=c)
            lc = lasso(data, cols=c)
            print("NN SCORE: {}, KRR SCORE: {}, LASSO SCORE: {}".format(avg_nnc, krrc, lc))

            print("\nWITHOUT APPROXIMATION")
            avg_nnnc = 0
            for i in range(10):
                avg_nnnc += 100
                #avg_nnnc += nn(data, cols=no_crude_c)/10
            krrnc = kernel_regression(data, cols=no_crude_c)
            lnc = lasso(data, cols=no_crude_c)
            print("NN SCORE: {}, KRR SCORE: {}, LASSO SCORE: {}\n\n".format(avg_nnnc,
                            krrnc, lnc))

            nn_scores.append(avg_nnc)
            nc_nn_scores.append(avg_nnnc)
            krr_scores.append(krrc)
            nc_krr_scores.append(krrnc)
            lasso_scores.append(lc)
            nc_lasso_scores.append(lnc)

            if(avg_nnc < best_score):
                best_score = avg_nnc
                best_col_combination = c
                best_model = "NN"
            if(avg_nnnc < best_score):
                best_score = avg_nnc
                best_col_combination = no_crude_c
                best_model = "NN"

            if(krrc < best_score):
                best_score = krrc
                best_col_combination = c
                best_model = "KRR"
            if(krrnc < best_score):
                best_score = krrc
                best_col_combination = no_crude_c
                best_model = "KRR"

            if(lc < best_score):
                best_score = lc
                best_col_combination = c
                best_model = "LASSO"
            if(lnc < best_score):
                best_score = lc
                best_col_combination = no_crude_c
                best_model = "LASSO"
            
    print("BEST COLUMN COMBINATION: {}, BEST SCORE: {}, BEST MODEL: {}".format(
               best_col_combination, best_score, best_model))

    output_data = pd.DataFrame({
        'columns': col_combs,
        'lasso rmse': lasso_scores,
        #'nn scores': nn_scores,
        'krr rmse': krr_scores,
        'lasso rmse no crude': nc_lasso_scores,
        #'nn scores no crude': nc_nn_scores,
        'krr rmse no crude': nc_krr_scores,
    })
    output_data.to_csv("./model_fittings.csv")

    if(best_model == 'KRR'):
        print(best_col_combination)
        krr_model = KR(kernel='rbf')
        krr_model.fit(data[best_col_combination], data["formation_energy_per_atom"])
        preds = krr_model.predict(data[best_col_combination])

        fig, ax = plt.subplots()
        ax.scatter(data["num_cu"], preds, label="Predicted Energy", marker='s', edgecolor='k',
                   color='b')
        ax.scatter(data["num_cu"], data["formation_energy_per_atom"], label="Target Energy",
                   marker='h', edgecolor='k', color='y')
        ax.set(xlabel="Number of Copper Atoms in Unit Cell",
               ylabel="Formation Energy Per Atom (eV)",
               title="Predictions From Best Model")
        ax.legend(loc='best')
        plt.show()
        exit(1)

    
    ##print("WITH APPROXIMATION")
    ##nn(data, cols=["crude", "a", "b", "c", "alpha", "beta", "gamma", "density"])
    ##kernel_regression(data, cols=["crude", "a", "b", "c", "alpha", "beta", "gamma", "density"])
    ##lasso(data, cols=["crude", "a", "b", "c", "alpha", "beta", "gamma", "density"])

    ##print("\nWITHOUT APPROXIMATION")
    ##nn(data, cols=["a", "b", "c", "alpha", "beta", "gamma", "density"])
    ##kernel_regression(data, cols=["a", "b", "c", "alpha", "beta", "gamma", "density"])
    ##lasso(data, cols=["a", "b", "c", "alpha", "beta", "gamma", "density"])

    #num = len(data["pretty_formula"].values)
    #fig, ax = plt.subplots()
    ##ax.scatter(data["num_cu"], data["formation_energy_per_atom"])
    ##ax.scatter(data["num_cu"], data["crude"])
    #ax.scatter(data["num_sn"], data["formation_energy_per_atom"])

    #plt.show()


