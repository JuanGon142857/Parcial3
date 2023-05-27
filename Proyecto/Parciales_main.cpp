#include <cmath>
#include "Parciales.h" 

int main(){
    double xL = -30.;
    double xR = 70.;
    double beta = -1;
    double gamma = 2./3;
    double T = 25;
    int N = 2001;
    int m = 2500;
    double TOL = 1e-12;

    Parcial myParcial(-30, 70, -1, 2./3, 50, 2001, 5000, 1e-12);    
    //xL = limite izquierdo del intervalo espacial
    //xR = limite derecho del intervalor espacial
    //
    myParcial.Set_X();
    myParcial.CrankNicolson( "One Soliton.csv" );
    
}