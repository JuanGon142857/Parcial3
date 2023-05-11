#include "No_Lineal.h" 

int main() {
    No_Lineal myNOLineal(0, 120., 0, 0, 19, 50, 1e-5);
    //myNOLineal.Set_N( 19);
    myNOLineal.Set_x();
    myNOLineal.Solve();
    myNOLineal.Solution2csv("Ejercicio 6 libro.csv");

}