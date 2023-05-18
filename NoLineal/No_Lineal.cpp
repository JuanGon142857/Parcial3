#include <iostream>
#include <cmath>
#include <iomanip>
#include <fstream>

using namespace std;

#include "No_Lineal.h"

No_Lineal :: No_Lineal(double x0, double xf, double y0, double yf, int N, int MaxIter, double TOL) 
: Lineal(x0, xf, y0, yf, N){
    Set_MaxIter( MaxIter );
    Set_TOL( TOL );
};

void No_Lineal :: Set_MaxIter( int k ){
    if (k <= 0){
        cout << "El maximo de iteraciones ingresado no es valido, se redefinio como 1\n";
        MaxIter = 1;
    }
    else{
        MaxIter = k;
    }
}

void No_Lineal :: Set_TOL( double tol ){
    if (tol <= 0){
        cout << "El valor de tolerancia ingresa no es valido, se redefinio como 1E-8\n";
        TOL = 1E-8;
    }
    else{
        TOL = tol;
    }
}

int No_Lineal :: Get_MaxIter() const{
    return MaxIter;
}

double No_Lineal :: Get_TOL() const{
    return TOL;
}

const double No_Lineal :: f( double x, double y, double yp){
    double S = 1000.;
    double I = 625.;
    double E = 3.0e7;
    double l = 120.;
    double q = 100./12; //constantes que se usan en el ejercicio 6

    return (S / (E * I) * y + q * x * (x - l) / (2 * E * I)) / pow((1 + yp * yp), 3./2) ;
}

void No_Lineal :: Solve(){
    

    omega_i.clear();
    for (int i = 0; i <= N + 1; i++){
        omega_i.push_back(0);
    }

    //Paso 1
    double h = (Get_xf() - Get_x0()) / (Get_N() + 1); //Divide el intervalo en secciones
    omega_i[0] = y0;   omega_i[N + 1] = yf;

    //Paso 2
    for (size_t i = 1; i <= N; i++){
        omega_i[i] = y0 + i * (yf - y0) / (xf - x0) * h;
    }

    std :: vector <double> a_i(N + 2);
    std :: vector <double> b_i(N + 2);
    std :: vector <double> c_i(N + 2);
    std :: vector <double> d_i(N + 2);

    std :: vector <double> l_i(N + 2);
    std :: vector <double> u_i(N + 2);
    std :: vector <double> v_i(N + 2);
    std :: vector <double> z_i(N + 2);

    //Paso 3
    size_t k = 1;

    bool Completado = false;

    //Paso 4    Pasos del 5 al 16 aplican el metodo de Newton 
    while ((!Completado) && (k <= MaxIter)){

        //Paso 5
        double x = x0 + h;
        double t = (omega_i[2] - y0) / (2 * h);

        a_i[1] = 2 + h * h * (f(x, omega_i[1] + h, t) - f(x, omega_i[1], t)) / h; 
        b_i[1] = -1 + (h / 2) * (f(x, omega_i[1], t + h) - f(x, omega_i[1], t)) / h; 
        d_i[1] = -(2 * omega_i[1] - omega_i[2] - y0 + h * h * f(x, omega_i[1], t));

        //Paso 6
        for (size_t i = 2; i <= N - 1; i++){
            x = x0 + i * h;
            t = (omega_i[i + 1] - omega_i[i - 1]) / (2 * h);
            a_i[i] = 2 + h * h * (f(x, omega_i[i] + h, t) - f(x, omega_i[i], t)) / h;
            b_i[i] = -1 + (h / 2) * (f(x, omega_i[i], t + h) - f(x, omega_i[i], t)) / h; 
            c_i[i] = -1 - (h / 2) * (f(x, omega_i[i], t + h) - f(x, omega_i[i], t)) / h; 
            d_i[i] = -(2 * omega_i[i] - omega_i[i + 1] - omega_i[i - 1] + h * h * f(x, omega_i[i], t));
        }

        //Paso 7
        x = xf - h;
        t = (yf - omega_i[N - 1]) / (2 * h);
        a_i[N] = 2 + h * h * (f(x, omega_i[N] + h, t) - f(x, omega_i[N], t)) / h;
        c_i[N] = -1 - (h / 2) * (f(x, omega_i[N], t + h) - f(x, omega_i[N], t)) / h;
        d_i[N] = -(2 * omega_i[N] - omega_i[N - 1] - yf + h * h * f(x, omega_i[N], t));

        //Paso 8. Pasos del 8 al 12 resuelven la matriz tridiagonal mediante factorazacion LU.
        l_i[1] = a_i[1];
        u_i[1] = b_i[1] / a_i[1];
        z_i[1] = d_i[1] / l_i[1];

        //Paso 9
        for (size_t i = 2; i <= N - 1; i++){
            l_i[i] = a_i[i] - c_i[i] * u_i[i - 1];
            u_i[i] = b_i[i] / l_i[i];
            z_i[i] = (d_i[i] - c_i[i] * z_i[i - 1]) / l_i[i];
        }

        //Paso 10
        l_i[N] = a_i[N] - c_i[N] * u_i[N - 1];
        z_i[N] = (d_i[N] - c_i[N] * z_i[N - 1]) / l_i[N];

        //Paso 12. Aplica la correcion iterativa
        v_i[N] = z_i[N];
        omega_i[N] = omega_i[N] + v_i[N];

        //Paso 13. Pasos 13 y 14 determinan si se ha alcanzado convergencia
        for (size_t i = N - 1; i >= 1; i--){
            v_i[i] = z_i[i] - u_i[i] * v_i[i + 1];
            omega_i[i] = omega_i[i] + v_i[i];
        }

        double suma_v = 0;
        for (auto& n : v_i){
            suma_v += n * n;
        }
        
        if (suma_v <= TOL * TOL){
            bool Completado = true;
        }
        k++;
    }
    if (!Completado){
    std :: cout << "Se excedio el maximo numero de iteraciones \n";
    }
    /*
    for (size_t i = 0; i <= N + 1; i++){
        omega_i[i] = cos(x_i[i]);
    }*/
}