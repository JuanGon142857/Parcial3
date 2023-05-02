#include <iostream>

#include "DiferenciasFinitas.h" 

int main(){
    DiferenciasFinitas myDiferenciasFinitas;
    myDiferenciasFinitas.Set_x0( float(1) );
    myDiferenciasFinitas.Set_xf( float(2) );
    myDiferenciasFinitas.Set_y0( float(1) );
    myDiferenciasFinitas.Set_yf( float(2) );
    myDiferenciasFinitas.Set_N( 9 );

    //myDiferenciasFinitas.Set_x();

    std :: vector <float> vi = myDiferenciasFinitas.Get_x();
    std :: cout << "Aca" << "\n ";
    myDiferenciasFinitas.Solve_Lineal();
    myDiferenciasFinitas.Imprimir_omega();
    //std :: cout << "Hi";

    return 0; 
} 
/*
int main(){
    
    int N = 9;
    float a = 1;
    float b = 2;
    float alpha = 1;
    float beta = 2;
    std :: vector <float> omega = Lineal(a, b, alpha, beta, p, q, r, N);
    for (size_t i = 0; i < N + 2; i++){
    std :: cout << omega[i] << "\n";
    }
    int N = 19;
    float a = 1;
    float b = 3;
    float alpha = 17;
    float beta = float(43) / 3;
    float TOL = 1E-8;
    int Max_iter = 10;
    std :: vector <float> omega = No_Lineal(a, b, alpha, beta, f, N, TOL, Max_iter);
    for (size_t i = 0; i < N + 2; i++){
    std :: cout << omega[i] << "\n";
    }

    return 0;
}*/