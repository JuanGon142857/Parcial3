#include <iostream>
#include <cmath>
#include <iomanip>
#include <fstream>

#include "Lineal.h" 

using namespace std;

Lineal :: Lineal(double x0, double xf, double y0, double yf, int N){
    Set_x0(x0);
    Set_xf(xf);
    Set_y0(y0);
    Set_yf(yf);
    Set_N(N);
}

void Lineal :: Set_x0(double a){
    x0 = a;
}

void Lineal :: Set_xf(double b){
    if (b == x0){
        std :: cout << "Se tiene xf = x0, se redefinio xf = x0 + 1 para evitar esto\n";
        xf = b + 1;
    }
    else{
        xf = b;
    }
    
}

void Lineal :: Set_y0(double alpha){
    y0 = alpha;
}

void Lineal :: Set_yf(double beta){
    yf = beta;
}

void Lineal :: Set_N(int n){
    if (n >= 2){
        N = n;
    }
    else{
        std :: cout << "El valor ingresado es demasiado pequeño, n se redefenio como 2 \n";
        N = 2;
    }
}

double Lineal :: Get_x0() const{
    return x0;
}

double Lineal :: Get_xf() const{
    return xf;
}

double Lineal :: Get_y0() const{
    return y0;
}

double Lineal :: Get_yf() const{
    return yf;
}

int Lineal :: Get_N() const{
    return N;
}


vector <double> Lineal :: Get_omega() const{
    return omega_i;
}

vector <double> Lineal :: Get_x() const{
    return x_i;
};

const double Lineal :: p(double x) {
    double delta = 0.3;
    return 0;
}

const double Lineal :: q(double x){
    double S = 1000.;
    double I = 625.;
    double E = 3e7;
    return S / (E * I);
}

const double Lineal :: r(double x){
    double S = 1000.;
    double I = 625.;
    double E = 3e7;
    double q = 100. / 12;
    double l = 120.;

    return q * x / (2 * E * I) * (x - l);
}

void Lineal :: Set_x() {
    x_i.clear();
    for (int i = 0; i <= N + 1; i++){
        x_i.push_back(x0 + (xf - x0) / (N + 1) * i);
    }
};


void Lineal :: Solve(){
    omega_i.clear();
    for (int i = 0; i <= N + 1; i++){
        omega_i.push_back(0);
    }
    
    //Paso 1
    double h = (xf - x0) / (N + 1); //Divide el intervalo en secciones
    double x = x0 + h;
    
    std :: vector <double> a_i(N + 2);   a_i[1] = 2 + h * h * q(x);
    std :: vector <double> b_i(N + 2);   b_i[1] = -1 + h / 2 * p(x);
    std :: vector <double> c_i(N + 2);
    std :: vector <double> d_i(N + 2);   d_i[1] = -h * h * r(x) + (1 + (h / 2) * p(x)) * y0;

    //Paso 2
    for (size_t i = 2; i <= N - 1; i++){
        x = x0 + i * h;
        a_i[i] = 2 + h * h * q(x);
        b_i[i] = -1 + (h / 2) * p(x);
        c_i[i] = -1 - (h / 2) * p(x);
        d_i[i] = -h * h * r(x);
    }

    //Paso 3
    x = xf - h;
    a_i[N] = 2 + h * h * q(x);
    c_i[N] = -1 - (h / 2) * p(x);
    d_i[N] = -h * h * r(x) + (1 - (h / 2) * p(x)) * yf;

    //Paso 4. Pasos del al 8 son factorizacion LU del sistema tridiagonal
    std :: vector <float> l_i(N + 2);   l_i[1] = a_i[1];
    std :: vector <float> u_i(N + 2);   u_i[1] = b_i[1] / a_i[1];
    std :: vector <float> z_i(N + 2);   z_i[1] = d_i[1] / l_i[1];

    //Paso 5
    for (size_t i = 2; i <= N - 1; i++){
        l_i[i] = a_i[i] - c_i[i] * u_i[i - 1];
        u_i[i] = b_i[i] / l_i[i];
        z_i[i] = (d_i[i] - c_i[i] * z_i[i - 1]) / l_i[i];
    }

    //Paso 6
    l_i[N] = a_i[N] - c_i[N] * u_i[N - 1];
    z_i[N] = (d_i[N] - c_i[N] * z_i[N - 1]) / l_i[N];

    //Paso 7
    omega_i[0] = y0;   omega_i[N + 1] = yf;    omega_i[N] = z_i[N];

    //Paso 8
    for (size_t i = N - 1; i >= 1; i--){
        omega_i[i] = z_i[i] - u_i[i] * omega_i[i + 1];
    }     
}

void Lineal :: Solution2csv( const char* filename ){
    ofstream Solucion( filename );
    if (!Solucion){
        cerr << "No se pudo crear el archivo \n";
        exit( 1 );
    }
    else{
        Solucion << "x,y\n";
        for (int i = 0; i <= Get_N() + 1; i++){
            Solucion << setw(6) << setprecision(5) << fixed << x_i[i] << "," << setw(8) << setprecision(5) << fixed << omega_i[i] << "\n"; 
        }
        Solucion.close();
    }
}