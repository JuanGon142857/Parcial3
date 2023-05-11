#include "No_Lineal.h" 

int main() {
    No_Lineal myNOLineal(0, 120., 0, 0, 19, 50, 1e-5); //x0, xf, y0, yf, N - 1 intervalos, numero maximo de iteraciones del metodo iterativo, tolerancia para definir si se alcazò la convergencia
    //myNOLineal.Set_N( 19);
    myNOLineal.Set_x(); //Encuentra los puntos x
    myNOLineal.Solve(); //Encuentra las aproximasciones y
    myNOLineal.Solution2csv("Ejercicio 6 libro.csv");

}