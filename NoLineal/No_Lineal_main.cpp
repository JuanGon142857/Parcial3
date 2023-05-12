#include <cmath>
#include "No_Lineal.h" 

int main() {
    No_Lineal myNOLineal(150, 150 + 4.48 * 1, 1.26 , 1.26, 10000, 1000, 1e-4); //x0, xf, y0, yf, N - 1 intervalos, numero maximo de iteraciones del metodo iterativo, tolerancia para definir si se alcazò la convergencia
    //myNOLineal.Set_N( 19);
    myNOLineal.Set_x(); //Encuentra los puntos x
    myNOLineal.Solve(); //Encuentra las aproximasciones y
    myNOLineal.Solution2csv("Duffin oscilator.csv");

}