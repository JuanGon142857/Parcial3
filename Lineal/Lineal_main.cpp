#include <cmath>

#include "Lineal.h" 

int main() {
    Lineal myLineal(0, 8.0 * 3.141592 / 1, 0, 0, 999); //x0, xf, y0, yf, N - 1 intervalos
    myLineal.Set_x();   //Encuentra los puntos x
    myLineal.Solve();   //Encuentra las aproximasciones y
    myLineal.Solution2csv("Oscilador forzado sin amoritugamiento resonancia.csv"); //El argumento es el nombre del archivo donde se guardaran los resultados

}