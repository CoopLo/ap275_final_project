from matplotlib import pyplot as plt
import pandas as pd

if __name__ == '__main__':
    data = pd.read_csv("./model_fittings.csv")
    #print(data["total_magnetization"])
    print(len(data))
    fig, ax = plt.subplots()
    ax.scatter(range(len(data)), data["krr rmse"], label="Crude Approximation", marker='s',
               edgecolor='k', color='r')
    ax.scatter(range(len(data)), data["krr rmse no crude"], label="No Crude Approximation",
               edgecolor='k', marker='p', color='g')
    ax.legend(loc='best')
    ax.set(xlabel="Combinations of Input Features", ylabel="Average MSE of Predictions",
           title="Comparison of KRR MSE Over All Combinations of Features")
    ax.set(xticks=[])
    plt.show()
    pass

