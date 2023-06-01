#include <cmath>
#include <complex>
#include "Parciales.h" 

using namespace std;
//Función de condiciones iniciales

complex <double> f( double);

double xL = -30.;
double xR = 70.;
double beta0 = -1;
double gamnma = 2./3;
double T = 25;
int N = 2001;
int m = 2500;
double TOL = 1e-12;

int main(){
    Parcial myParcial(xL, xR, beta0, gamnma, T, N, m, TOL, f);    
    //xL = limite izquierdo del intervalo espacial
    //xR = limite derecho del intervalor espacial
    //
    myParcial.Set_X();
    myParcial.CrankNicolson( "One Soliton.csv" );
    
}

complex <double> f( double x) {
    double a = 1;
    double c = 1;
    return sqrt( 2 * a / gamnma ) * 1. / cosh(sqrt(-2 * a / beta0) * x) * exp( 1i * (-c / beta0) * x);
};