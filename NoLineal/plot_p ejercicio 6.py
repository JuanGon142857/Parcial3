import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import csv

data = pd.read_csv("./Ejercicio 6 libro.csv")
dataL = pd.read_csv("./../Lineal/Ejercicio 7 libro.csv")

def f(x):
    c1 = 7.7042537e4
    c2 = 7.9207462e4
    a = 2.3094010e-4
    b = -4.1666666e-3
    c = -1.5625e5
    l = 120
    return c1 * np.exp(a * x) + c2 * np.exp(-a * x) + b * (x - l) * x + c 

plt.plot(data.x.values, (f(data.x.values) - f(0))*1000, c = 'k')
plt.scatter(data.x.values, dataL.y.values*1000, c = 'b', s = 500, marker= "*")
plt.scatter(data.x.values, data.y.values*1000, c = 'r', s = 100)

plt.xticks(np.linspace(0, 120, 5), fontsize = 30)
plt.yticks(fontsize = 30)
plt.legend(["Numérica", "Método lineal", "Método no lineal"], fontsize = 30)
plt.xlabel("posición (pulgadas)", fontsize = 30)
plt.ylabel("deflexión (pulgadas / 1000)", fontsize = 30)
plt.show()
print(data.y.values)
print(f(data.x.values))
print(abs(data.y.values - f(data.x.values)))

