#include <cmath>
#include "Parciales.h" 

int main(){
    Parcial myParcial(-30, 70, -1, 2./3, 25, 2001, 2500, 1e-12);    //l, T, m, N, alpha
    myParcial.Set_X();
    myParcial.CrankNicolson( "One Soliton.csv" );
    
}