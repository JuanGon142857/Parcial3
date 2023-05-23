#include <cmath>
#include "Parciales.h" 

int main(){
    Parcial myParcial(-30, 70, -1, 2./3, 70, 2001, 5000, 1e-12);    //l, T, m, N, alpha
    myParcial.Set_X();
    myParcial.CrankNicolson( "One Soliton.csv" );
    
}