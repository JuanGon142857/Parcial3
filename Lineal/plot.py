import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import csv

data = pd.read_csv("./Oscilador forzado sin amoritugamiento resonancia.csv")

def f(x):
    return np.exp(-0.4/2 * x)

plt.plot(data.x.values, data.y.values)
#plt.plot(data.x.values, f(data.x.values), linestyle = 'dashed')
#plt.legend(["Numerico", "Amplitud"])
plt.show()
#print(data.x.values)

