import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import csv

data = pd.read_csv("./ejercicio 1 Libro.csv")

y = data.loc[0].values
y = y[1:]

x = np.linspace(0.5, 2., 4)

print(x)

def f(x,t):
    return np.exp(-(np.pi / 2) ** 2 * t) * np.sin(np.pi/2 * x) 

plt.plot(x,y)
plt.plot(x, f(x, 0.1), linestyle = 'dashed')
plt.show()
#plt.plot(data.x.values, f(data.x.values), linestyle = 'dashed')
#plt.legend(["Numerico", "Amplitud"])
#print(data.x.values)

