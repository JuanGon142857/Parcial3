#define _USE_MATH_DEFINES

#include <iostream>
#include <cmath>
#include <iomanip>
#include <fstream>

#include "Parciales.h" 

using namespace std;

Parcial :: Parcial(float l, float T, int m, int N, float alpha){
    Set_l(l);
    Set_T(T);
    Set_m(m);
    Set_N(N);
    Set_alpha(alpha);
}

void Parcial :: Set_l(float b){
    if (b <= 0){
        std :: cout << "l no puede ser no positivo, se redefinio l =  1 para evitar esto\n";
        l = 1;
    }
    else{
        l = b;
    }
}

void Parcial :: Set_T(float beta){
    if (beta <= 0){
        std :: cout << "T no puede ser no positivo, se redefinio l =  1 para evitar esto\n";
        T = 1;
    }
    else{
        T = beta;
    }
}

void Parcial :: Set_m(int M){
    if (M >= 3){
        m = M;
    }
    else{
        std :: cout << "El m ingresado es demasiado pequeño, m se redefenio como 3 \n";
        m = 3;
    }
}

void Parcial :: Set_N(int n){
    if (n >= 1){
        N = n;
    }
    else{
        std :: cout << "El N ingresado es demasiado pequeño, n se redefenio como 1 \n";
        N = 1;
    }
}

void Parcial :: Set_alpha(float a){
    alpha = a;
}


float Parcial :: Get_l() const{
    return l;
}

float Parcial :: Get_T() const{
    return T;
}

int Parcial :: Get_N() const{
    return N;
}

int Parcial :: Get_m() const{
    return m;
}

float Parcial :: Get_alpha() const{
    return alpha;
}

vector <float> Parcial :: Get_omega() const{
    return omega_i;
}

vector <float> Parcial :: Get_x() const{
    return x_i;
};

void Parcial :: Set_x() {
    x_i.clear();
    for (int i = 0; i <= m; i++){
        x_i.push_back(i * l / m);
    }
};

float Parcial :: f( float x) {
    return 1 / (M_PI * M_PI)  * sin(x * M_PI / l);
};

void Parcial :: CrankNicolson( const char* filename ){
    omega_i.clear();
    for (int i = 0; i <= m; i++){
        omega_i.push_back(0);
    }
    
    //Paso 1
    float h = l / m;    
    float k = T / N;
    float lambda = alpha * alpha * k / (h * h);
    omega_i[m] = 0;

    //Paso 2
    for (int i = 1; i <= m - 1; i++){
        omega_i[i] = f(i * h);
    }

    //Factorizacion LU del paso 3 al 11

    //Paso 3
    std :: vector <float> l_i(m + 1);   l_i[1] = 1. + lambda;
    std :: vector <float> u_i(m + 1);   u_i[1] = - lambda / (2 * l_i[1]);

    //Paso 4
    for (int i = 2; i <= m - 2; i++){
        l_i[i] = 1 + lambda + lambda * u_i[i - 1] / 2.;
        u_i[i] = -lambda / (2 * l_i[i]); 
    }

    //Paso 5
    l_i[m - 1] = 1 + lambda + lambda * u_i[m - 2] / 2.;

    std :: vector <float> z_i(m + 1);

    //Paso 6
    for (int j = 1; j <= N; j++){
        
        //Paso 7
        int t = j * k;
        z_i[1] = ((1 - lambda) * omega_i[1] + lambda / 2. * omega_i[2]) / l_i[1];

        //Paso 8
        for (int i = 2; i <= m - 1; i++){
            z_i[i] = ((1 - lambda) * omega_i[i] + lambda / 2. * (omega_i[i + 1] + omega_i[i - 1] + z_i[i - 1])) / l_i[i];
            u_i[i] = -lambda / (2 * l_i[i]); 
        }

        //Paso 9
        omega_i[m - 1] = z_i[m - 1];

        //Paso 10
        for (int i = m - 2; i >= 1; i--){
            omega_i[i] = z_i[i] - u_i[i] * omega_i[i + 1];
        }

        if (j == 1){
            ofstream Solucion( filename);
                if (!Solucion){
                cerr << "No se pudo crear el archivo \n";
                exit( 1 );
            }
            else{
                Solucion << "T";
                for (int i = 1; i <= Get_m(); i++){
                    Solucion << ",x" << i; 
                }
                Solucion << "\n";
                Solucion.close();
            }
        }

        //Devolver
        ofstream Solucion( filename, ios :: app);
        if (!Solucion){
            cerr << "No se pudo crear el archivo \n";
            exit( 1 );
        }
        
        else{
            Solucion << j * k;
            for (int i = 1; i <= Get_m(); i++){
                Solucion << "," << setw(9) << setprecision(8) << fixed << omega_i[i]; 
            }
            Solucion << "\n";
            Solucion.close();
        }

    }



}

void Parcial :: Solution2csv( const char* filename ){
    ofstream Solucion( filename, ios :: app);
    if (!Solucion){
        cerr << "No se pudo crear el archivo \n";
        exit( 1 );
    }
    
    else{
        for (int i = 1; i <= Get_m(); i++){
            Solucion << setw(6) << setprecision(5) << fixed << omega_i[i] << ","; 
        }
        Solucion << "\n";
        Solucion.close();
    }
}