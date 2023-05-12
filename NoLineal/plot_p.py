import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import csv

data = pd.read_csv("./Duffin oscilator.csv")

def f(x):
    return x* np.sin(x) / 2


plt.plot(data.x.values - 150, data.y.values)
#plt.plot(data.x.values, f(data.x.values), linestyle = 'dashed')
plt.legend(["Numerico", "Exacto"])
plt.show()
#print(data.x.values)

