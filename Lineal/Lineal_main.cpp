#include "Lineal.h" 

int main() {
    Lineal myLineal(0, 50, 0, 0, 49);
    myLineal.Set_x();
    myLineal.Solve();
    myLineal.Solution2csv("Ejercicio 8 libro.csv");

}