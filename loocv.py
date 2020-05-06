import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
from sklearn.neural_network import MLPRegressor
from sklearn.kernel_ridge import KernelRidge as KR
from sklearn.linear_model import Lasso
from sklearn.model_selection import LeaveOneOut, cross_val_score
from ast import literal_eval
from itertools import combinations
#from numba import jit
import time
#import autograd as ad
#from autograd import numpy as np


def nn():#data, cols=["num_cu", "num_sn", "crude"]):
    model = MLPRegressor(activation='relu', solver='lbfgs', learning_rate='adaptive',
                 hidden_layer_sizes=(20), verbose=False,
                 max_iter=20000, tol=1e-9)
    #results = leave_one_out(data, model, cols)
    #print("NN AVERAGE RMSE: {}, IN {} SECONDS".format(results[0], results[1]))
    #return results[0]
    return model


def kernel_regression():#data, cols=["crude"]):
    sigma = 10.
    def kernel(x1, x2):
        return np.exp(-np.linalg.norm(x1-x2)**2/(2*sigma**2))

    def cos_sim(x1, x2):
        return np.dot(x1, x2)/(np.linalg.norm(x1)*np.linalg.norm(x2))

    model = KR(kernel=cos_sim, alpha=0.0000000001)
    #return leave_one_out(data, model, cols)
    #results = leave_one_out(data, model, cols)
    #print("KRR AVERAGE RMSE: {}, IN {} SECONDS".format(results[0], results[1]))
    #return results[0]
    return model


def lasso():#data, cols=["crude"]):
    model = Lasso()
    #results = leave_one_out(data, model, cols)
    #print("LASSO AVERAGE RMSE: {}, IN {} SECONDS".format(results[0], results[1]))
    #return results[0]
    #return leave_one_out(data, model, cols)
    return model


def leave_one_out(model, data, cols):#data, model, cols):
    print(cols)
    target = 'formation_energy_per_atom'
    loo = LeaveOneOut()
    rmse = []
    cv_scores = []
    start = time.time()
    for train_index, test_index in loo.split(data[cols]):
        tmodel = lasso() if(model=='lasso') else kernel_regression() if(model=='kernel') else \
                 nn() if(model=='nn') else None

        # Split dataset
        x_train, x_test = data[cols].loc[train_index], data[cols].loc[test_index]
        y_train, y_test = data[target].loc[train_index], data[target].loc[test_index]

        # Fit model
        tmodel.fit(x_train, y_train)

        # Get metrics
        pred_y = tmodel.predict(x_test)
        rmse.append((pred_y[0] - y_test.values[0])**2)
        #cv_scores.append(cross_val_score(model, x_test, y_test, cv=1))

    #print("AVERAGE RMSE: {}, IN {} SECONDS".format(np.sqrt(np.average(rmse)),
    #                                                 time.time()-start))
    #print(np.average(cv_scores))
    #exit(1)
    print("{} AVERAGE RMSE: {}, TIME: {}".format(model, np.sqrt(np.average(rmse)),
                                                 time.time()-start))
    return np.sqrt(np.average(rmse))#, time.time() - start
    #print(y_test.values.reshape(-1,1), pred_y.reshape(-1,1))
    #return model.score(y_test.values.reshape(-1,1), pred_y.reshape(-1,1))


if __name__ == '__main__':
    data = pd.read_csv("./cleaned_bronze_calculations.csv")
    data["cu_frac"] = data["num_cu"]/data["num_atoms"]
    print(data.columns)
    cols = ['density', 'volume', #'total_magnetization',
            'num_cu', 'num_sn']#, 'lattice', 'angles']

    #for col in cols:
        #if(('Unnamed' in col) or (not isinstance(data[col].values[0], float)) or
        #    col in ['formation_energy_per_atom', 'elasticity']):
        #    print(col)
        #    continue
    #    print(col)
    #    fig, ax = plt.subplots()
    #    ax.scatter(data[col], data["formation_energy_per_atom"])
    #    ax.set(title=col)
    #    plt.savefig("./trend/{}.png".format(col))
    #    plt.close()
    #exit(1)
    #fig, ax = plt.subplots()
    #ax.scatter(data["cu_frac"], data["formation_energy_per_atom"])
    #ax.scatter(data["cu_frac"], data["crude"])
    #ax.scatter(data["cu_frac"], data["not_so_crude"])
    #plt.show()
    #exit(1)
    #print(data[["pretty_formula", "crude", "formation_energy_per_atom"]])
    #exit(1)
    data["diff"] = data["formation_energy_per_atom"] - data["crude"]
    #print(data.columns)
    #print(data[data["pretty_formula"]=="Cu"])
    #exit(1)

    cols = ['density', 'volume', #'total_magnetization',
            'num_cu', 'num_sn', 'lattice', 'angles']
    target = 'formation_energy_per_atom'

    best_score = 99999
    best_col_combination = None
    best_model = None

    nn_scores = []
    krr_scores = []
    lasso_scores = []
    nc_nn_scores = []
    nc_krr_scores = []
    nc_lasso_scores = []

    col_combs = []
    for i in range(1, len(cols)+1):
        if(i > 2):
            continue
        for c in list(combinations(cols, i)):
            #col_combs.append(c)
            c = ['density', 'volume', 'lattice', 'angles']

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
            c.append("not_so_crude")
            #c.append("crude")

            #print("COLUMNS: {}".format(c))
            print("COLUMNS: {}".format(no_crude_c))
            print("WITH APPROXIMATION")
            avg_nnc = 0
            for i in range(1):
                avg_nnc += leave_one_out('nn', data, cols=c)/1
                #avg_nnc += 100
            krrc = leave_one_out('kernel', data, cols=c)
            #krrc = 100
            #lc = leave_one_out('lasso', data, cols=c)
            lc = 100
            print("NN SCORE: {}, KRR SCORE: {}, LASSO SCORE: {}".format(avg_nnc, krrc, lc))

            print("\nWITHOUT APPROXIMATION")
            avg_nnnc = 0
            for i in range(1):
                avg_nnnc += leave_one_out('nn', data, cols=no_crude_c)/1
                #avg_nnnc += 100
            krrnc = leave_one_out('kernel', data, cols=no_crude_c)
            #krrnc = kernel_regression(data, cols=no_crude_c)
            #krrnc = 100
            #lnc = leave_one_out('lasso', data, cols=no_crude_c)
            #lnc = lasso(data, cols=no_crude_c)
            lnc = 100
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
            break
        break
            
    print("BEST COLUMN COMBINATION: {}, BEST SCORE: {}, BEST MODEL: {}".format(
               best_col_combination, best_score, best_model))

    #output_data = pd.DataFrame({
    #    'columns': col_combs,
    #    'lasso rmse': lasso_scores,
    #    'nn scores': nn_scores,
    #    'krr rmse': krr_scores,
    #    'lasso rmse no crude': nc_lasso_scores,
    #    'nn scores no crude': nc_nn_scores,
    #    'krr rmse no crude': nc_krr_scores,
    #})
    #output_data.to_csv("./model_fittings.csv")

    #if(best_model == 'KRR'):
    #    #print(best_col_combination)
    #    model = KR(kernel='rbf')
    #    #model.fit(data[best_col_combination], data["formation_energy_per_atom"])
    #    #preds = krr_model.predict(data[best_col_combination])
    #elif(best_model == 'LASSO'):
    #    model = Lasso()
    #elif(best_model == 'NN'):
    #    model = MLPRegressor(activation='relu', solver='lbfgs', learning_rate='adaptive',
    #                hidden_layer_sizes=(20), verbose=False,
    #                max_iter=20000, tol=1e-9)

    preds = []
    loo = LeaveOneOut()
    rmse = []
    if('not_so_crude' not in best_col_combination):
        bc = list(best_col_combination)
        bc.append('not_so_crude')
        best_col_combination = np.array(bc)
        #list(best_col_combination).append('crude')
    #if('crude' not in best_col_combination):
    #    bc = list(best_col_combination)
    #    bc.append('crude')
    #    best_col_combination = np.array(bc)
    #    #list(best_col_combination).append('crude')
    print(best_col_combination)

    for train_index, test_index in loo.split(data[best_col_combination]):
        cmodel = lasso() if(best_model=='LASSO') else kernel_regression() \
                         if(best_model=='KRR') else nn() if(best_model=='NN') else None

        # Split dataset
        x_train = data[best_col_combination].loc[train_index]
        x_test = data[best_col_combination].loc[test_index]

        y_train = data[target].loc[train_index]
        y_test = data[target].loc[test_index]

        # Fit model
        cmodel.fit(x_train, y_train)

        # Get metrics
        pred_y = cmodel.predict(x_test)
        rmse.append((pred_y[0] - y_test.values[0])**2)
        preds.append(pred_y)
    print("CRUDE RMSE: {}".format(np.sqrt(np.mean(rmse))))

    ncpreds = []
    loo = LeaveOneOut()
    nrmse = []
    print(best_col_combination[:-1])
    for train_index, test_index in loo.split(data[best_col_combination[:-1]]):
        nmodel = lasso() if(best_model=='LASSO') else kernel_regression() \
                         if(best_model=='KRR') else nn() if(best_model=='NN') else None

        # Split dataset
        x_train = data[best_col_combination[:-1]].loc[train_index]
        x_test = data[best_col_combination[:-1]].loc[test_index]

        y_train = data[target].loc[train_index]
        y_test = data[target].loc[test_index]

        # Fit model
        nmodel.fit(x_train, y_train)

        # Get metrics
        pred_y = nmodel.predict(x_test)
        nrmse.append((pred_y[0] - y_test.values[0])**2)
        ncpreds.append(pred_y)
    print("NO CRUDE RMSE: {}".format(np.sqrt(np.mean(nrmse))))


    fig, ax = plt.subplots()
    ax.scatter(data["cu_frac"], preds,
               label="Not So Crude Predicted Energy, RMSE: {0:.2e}".format(np.sqrt(np.mean(rmse))),
               marker='s', edgecolor='k', color='b')
    ax.scatter(data["cu_frac"], ncpreds,
               label="No Crude Estimate Energy, RMSE: {0:.2e}".format(np.sqrt(np.mean(nrmse))),
               marker='p', edgecolor='k', color='r')
    ax.scatter(data["cu_frac"], data["formation_energy_per_atom"], label="Target Energy",
               marker='h', edgecolor='k', color='y')
    ax.set(xlabel="Cu Fraction",
           ylabel="Formation Energy Per Atom (eV)",
       title="Predictions From Best Model Using Set Columns")#: {0}".format(best_col_combination))
    #ax.legend(loc='best')2
    ax.legend(loc='upper center')
    plt.show()
    #exit(1)

    
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


