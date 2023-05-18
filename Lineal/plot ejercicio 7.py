import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import csv

data = pd.read_csv("./Ejercicio 7 libro.csv")

def f(x):
    c1 = 7.7042537e4
    c2 = 7.9207462e4
    a = 2.3094010e-4
    b = -4.1666666e-3
    c = -1.5625e5
    l = 120
    return c1 * np.exp(a * x) + c2 * np.exp(-a * x) + b * (x - l) * x + c 

plt.subplot(1, 2, 1)
plt.scatter(data.x.values, data.y.values * 1000, c = 'r')
plt.plot(data.x.values, (f(data.x.values) - f(0)) * 1000)
plt.legend(["Numérico", "Exacto"], fontsize = 20)
plt.xlabel("posición (pulgadas)", fontsize = 20)
plt.ylabel("deflexión (pulgadas / 1000)", fontsize = 20)
plt.xticks( fontsize = 20)
plt.yticks(fontsize = 20)
plt.subplot(1, 2, 2)
plt.plot(data.x.values, abs(data.y.values - f(data.x.values) + f(0)))
plt.xlabel("posición (pulgadas)", fontsize = 20)
plt.ylabel("Error absoluto", fontsize = 20)
plt.xticks( fontsize = 20)
plt.yticks(fontsize = 20)
plt.show()
print(data.y.values)
print(f(data.x.values))
print(abs(data.y.values - f(data.x.values)))

