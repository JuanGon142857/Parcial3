import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import csv

data = pd.read_csv("./Ejercicio 1 libro.csv")

def f(x):
    return np.exp(2) / (np.exp(4) - 1) * (np.exp(2 * x) - np.exp(-2 * x)) + x

plt.scatter(data.x.values, f(data.x.values), c = 'r')
plt.plot(data.x.values, data.y.values)
plt.legend(["Numérico", "Exacto"], fontsize = 15)
plt.xticks([0, 0.25, 0.5, 0.75, 1], fontsize = 15)
plt.yticks(np.linspace(0, 2, 5), fontsize = 15)
plt.show()
print(data.y.values)
print(f(data.x.values))
print(abs(data.y.values - f(data.x.values)))