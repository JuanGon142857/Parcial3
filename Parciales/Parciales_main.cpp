#include <cmath>
#include "Parciales.h" 

int main(){
    Parcial myParcial(1., 1., 100, 100, 1./ 3.14159);    //l, T, m, N, alpha
    myParcial.Set_x();
    myParcial.CrankNicolson( "Heat equation.csv" );
    
}