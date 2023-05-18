import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import csv

data = pd.read_csv("./Ejercicio 5 libro.csv")

def f(x):
    return np.exp(-10 * x)

plt.scatter(data.x.values, data.y.values, c = 'r')
plt.plot(data.x.values, f(data.x.values))
plt.legend(["Numérico", "Exacto"], fontsize = 15)
plt.xticks(np.linspace(0, 1, 5), fontsize = 15)
plt.yticks(np.linspace(0, 1, 5), fontsize = 15)
plt.show()
print(data.y.values)
print(f(data.x.values))
print(abs(data.y.values - f(data.x.values)))

