import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import csv

data = pd.read_csv("./Ejercicio 4d libro.csv")

def f(x):
    return x ** 3 *  np.log(x)

plt.plot(data.x.values, f(data.x.values))
plt.scatter(data.x.values, data.y.values, c = 'r')

plt.xticks(np.linspace(1, 2, 5), fontsize = 20)
plt.yticks(np.round(np.linspace(0, np.log(256), 5),2), fontsize = 20)
plt.legend(["Exacto", "Numérico"], fontsize = 20)
plt.show()
print(np.round(data.y.values, 5))
print(np.round(f(data.x.values), 5))
print(np.round(abs(data.y.values - f(data.x.values)), 5 ))

