#include <cmath>

#include "Lineal.h" 

int main() {
    Lineal myLineal(0., 120., 0., 0., 19); //x0, xf, y0, yf, N - 1 intervalos
    myLineal.Set_x();   //Encuentra los puntos x
    myLineal.Solve();   //Encuentra las aproximasciones y
    myLineal.Solution2csv("Ejercicio 7 libro.csv"); //El argumento es el nombre del archivo donde se guardaran los resultados

}