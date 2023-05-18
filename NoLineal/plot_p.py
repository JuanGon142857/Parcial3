import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import csv

data = pd.read_csv("./Ejercicio 2 libro.csv")

def f(x):
    return 1 / (x + 3)

plt.plot(data.x.values, f(data.x.values))
plt.scatter(data.x.values, data.y.values, c = 'r')

plt.xticks(np.linspace(-1, 0, 5), fontsize = 20)
plt.yticks(np.linspace(0.3, 0.5, 5), fontsize = 20)
plt.legend(["Exacto", "Numérico"], fontsize = 20)
plt.show()
print(data.y.values)
print(f(data.x.values))
print(abs(data.y.values - f(data.x.values)))

