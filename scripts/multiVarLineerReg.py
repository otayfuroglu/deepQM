#
import numpy as np
import pandas as pd
from sklearn import linear_model


def genUniqRandomInt(rng, size):

    cond = True
    while cond:
        rnd = np.random.randint(rng, size=size)
        if len(np.unique(rnd)) == len(rnd):
            cond = False
    return rnd


def calcRMSE(predictions, targets):
    return np.sqrt(((predictions - targets) ** 2).mean())


df = pd.read_csv("Bootsrtrap_ANID3LIE.csv")
#  print(genUniqRandomInt(len(df), size=6))
#  rnd = np.random.randint(len(df), size=6)
#  print(len(np.unique(rnd)) == len(rnd))

labels = ["Beta", "Alpha", "Gamma", "RMSE_train", "RMSE_test"]
df_result = pd.DataFrame(columns=labels)

b_values = []
a_values = []
g_values = []
rmse_train_values = []
rmse_test_values = []

for _ in range(100):
    rnd = genUniqRandomInt(len(df), size=6)
    msk = np.array([item in rnd for item in range(len(df))])
    #  print(msk)
    #
    df_train = df[~msk]
    df_test = df[msk]
    #
    #  print(len(df_train))
    #  print(len(df_test))
    #
    X_train = df_train[["coul", "lj"]]
    y_train = df_train["Experimental"]
    #
    X_test = df_test[["coul", "lj"]]
    y_test = df_test["Experimental"]
    regr = linear_model.LinearRegression()
    regr.fit(X_train, y_train)
    #
    predictedLIE_train = regr.predict(X_train)
    predictedLIE_test = regr.predict(X_test)

    b, a = regr.coef_
    g = regr.intercept_

    rmse_train = calcRMSE(predictedLIE_train, y_train.to_numpy())
    rmse_test = calcRMSE(predictedLIE_test, y_test.to_numpy())

    b_values += [b]
    a_values += [a]
    g_values += [g]
    rmse_train_values += [ rmse_train ]
    rmse_test_values += [ rmse_test ]


df_result[labels[0]] = b_values
df_result[labels[1]] = a_values
df_result[labels[2]] = g_values
df_result[labels[3]] = rmse_train_values
df_result[labels[4]] = rmse_test_values



df_result.to_csv("reslut_Bootsrtrap_ANID3LIE.csv")
