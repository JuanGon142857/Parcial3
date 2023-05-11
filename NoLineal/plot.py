import pandas as pd
import matplotlib.pyplot as plt
import csv

data = pd.read_csv("./Ejercicio 6 libro.csv")

plt.plot(data.x.values, data.y.values)
#plt.plot(data.x.values, 1 / (data.x.values + 3))
#plt.legend(["Numerico", "Exacto"])
plt.show()
#print(data.x.values)

