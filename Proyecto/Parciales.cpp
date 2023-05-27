#define _USE_MATH_DEFINES

#include <iostream>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <complex>
#include <algorithm>

#include "Parciales.h" 

using namespace std;

Parcial :: Parcial(double xL, double xR, double Beta, double Gamma, double T, int N, int m, double TOL){
    Set_xL(xL);
    Set_xR(xR);
    Set_beta(Beta);
    Set_gamma(Gamma);
    Set_T(T);
    Set_N(N);
    Set_m(m);
    Set_TOL(TOL);
}

void Parcial :: Set_xL(double b){
    xL = b;
}

void Parcial :: Set_xR(double b){
    xR = b;
}

void Parcial :: Set_beta(double b){
    beta = b;
}

void Parcial :: Set_gamma(double b){
    gamma = b;
}

void Parcial :: Set_T(double beta){
    if (beta <= 0){
        std :: cout << "T no puede ser no positivo, se redefinio l =  1 para evitar esto\n";
        T = 1;
    }
    else{
        T = beta;
    }
}

void Parcial :: Set_N(int n){
    if (n >= 2){
        N = n;
    }
    else{
        std :: cout << "El N ingresado es demasiado pequeño, n se redefenio como 1 \n";
        N = 2;
    }
    h = (xR - xL) / (N - 1);
}

void Parcial :: Set_m(int M){
    if (M >= 2){
        m = M;
    }
    else{
        std :: cout << "El m ingresado es demasiado pequeño, m se redefenio como 3 \n";
        m = 2;
    }
    k = T / (m - 1);
}

void Parcial :: Set_TOL(double tol){
    TOL = tol;
}

double Parcial :: Get_T() const{
    return T;
}

int Parcial :: Get_N() const{
    return N;
}

int Parcial :: Get_m() const{
    return m;
}

vector <double> Parcial :: Get_X() const{
    return X;
};

void Parcial :: Set_X() {
    X.clear();
    for (int i = 0; i <= m; i++){
        X.push_back(i * h);
    }
};

//Funcion para sacar el inverso de matriz 2x2
vector <vector <double>> Parcial :: inv_mat(vector <vector <double>> A){
    double det = A[0][0] * A[1][1] - A[1][0] * A[0][1];
    if (A.size() == 2 && A[0].size() == 2){
        if (det != 0){
            vector <vector <double>> inverse(2,vector<double> (2));
            inverse[0][0] = 1 / det * A[1][1];
            inverse[1][1] = 1 / det * A[0][0];
            inverse[0][1] = -1 / det * A[0][1];
            inverse[1][0] = -1 / det * A[1][0];
            return inverse;
        }
        else{
            cout << "Esta matriz no es invertible \n";
        }
    }
    else{
        cout << "Se ingresó una matriz que no es 2x2 \n";
    }
    return vector<vector <double>>{0};
    
}; 

//Funcion para multiplicar 2 matrices 2x2
vector <vector <double>> Parcial :: mat_mat_mul(vector <vector <double>> A,vector <vector <double>> B){
    if (A[0].size() == B.size()){
        vector <vector <double>> C(A.size(), vector<double> (B[0].size()));
        for (int i = 0; i < B[0].size(); i++){
            for (int j = 0; j < A.size(); j++){
                double suma = 0;
                for (int k = 0; k < A[0].size(); k++){
                    suma = suma + A[i][k] * B[k][j];
                } 
                C[i][j] = suma;
            }
        }
        return C;
    }
    else{
        cout << "Estas matrices no pueden multiplicarse entre s] \n";
    }
    return vector<vector <double>>{0};
}

//Funcion para aplicar matriz 2x2 sobre vector 2x1
vector <double> Parcial :: mat_vec_mul(vector <vector <double>> M, vector <double> v){
    if (M[0].size() == v.size()){
        vector <double> vp(M.size());
        for (int i = 0; i < v.size(); i++){
            for (int j = 0; j < M.size(); j++){
                double suma = 0;
                for (int k = 0; k < M[0].size(); k++){
                    suma = suma + M[i][k] * v[k];
                } 
                vp[i] = suma;
            }
        }
        return vp;
    }
    else{
        cout << "Estas matrices no pueden multiplicarse entre sì \n";
    }
    return vector<double>{0};
}

//Funcion para sumar dos matrices
vector<vector <double>> Parcial :: mat_mat_sum(vector <vector <double>> A, vector <vector <double>> B){
    if ((A.size() == B.size()) && (A[0].size() == B[0].size())){
        vector <vector <double>> C = A;
        for (int i = 0; i < A.size(); i++){
            for (int j = 0; j < A[0].size(); j++){
                C[i][j] = A[i][j] + B[i][j];
            }
        }
        return C;
    }
    else{
        cout << "Estas matrices no pueden sumarse entre sì \n";
    }
    return vector<vector <double>>{0};
}

//Funcion para multiplicar matriz por una constante
vector<vector <double>> Parcial :: mat_const_mul(vector <vector <double>> A, double c){
    vector <vector <double>> cA = A;
    for (int i = 0; i < A.size(); i++){
        for (int j = 0; j < A[0].size(); j++){
            cA[i][j] = A[i][j] * c;
        }
    }
    return cA;
}

//Funcion para sumar 2 vectores
vector <double> Parcial :: vec_vec_sum(vector <double> A, vector <double> B){
    if ((A.size() == B.size())){
        vector <double> C = A;
        for (int i = 0; i < A.size(); i++){
            C[i] = A[i] + B[i];
        }
        return C;
    }
    else{
        cout << "Estas matrices no pueden sumarse entre sì \n";
    }
    return vector <double>{0};
}

//Funcion para multiplicar vector por una constante
vector <double> Parcial :: vec_const_mul(vector <double> A, double c){
    vector <double> cA = A;
    for (int i = 0; i < A.size(); i++){
        cA[i] = A[i] * c;
    }
    return cA;
}

//Función de condiciones iniciales
complex <double> Parcial :: f( double x) {
    double a = 10;
    double c = 1;
    return sqrt( 2 * a / gamma ) * 1. / cosh(sqrt(-2 * a / beta) * x) * exp( 1i * (-c / beta) * x);
};

//Funcion que define el sistema de ecuaciones no lineales
vector <double> Parcial :: F(int i){
    double modulo = pow((temp[i][0] + omega[i][0]), 2.) + pow((temp[i][1] + omega[i][1]), 2.);
    double REAL;
    double IMAG;
    if (i == 0){
        REAL = temp[i][0] - omega[i][0] + beta * k / (2. * h * h) * (2. * temp[i][1] + 2. * omega[i][1] - 2. * temp[i + 1][1] - 2. * omega[i + 1][1]) / 2. + k * gamma / 4. * modulo * (temp[i][1] + omega[i][1]) / 2.;
        IMAG = temp[i][1] - omega[i][1] - beta * k / (2. * h * h) * (2. * temp[i][0] + 2. * omega[i][0] - 2. * temp[i + 1][0] - 2. * omega[i + 1][0]) / 2. - k * gamma / 4. * modulo * (temp[i][0] + omega[i][0]) / 2.;
    }
    if (i > 0 && i < N - 1){
        REAL = temp[i][0] - omega[i][0] + beta * k / (2. * h * h) * (- temp[i - 1][1] - omega[i - 1][1] + 2. * temp[i][1] + 2. * omega[i][1] - temp[i + 1][1] - omega[i + 1][1]) / 2. + k * gamma / 4. * modulo * (temp[i][1] + omega[i][1]) / 2.;
        IMAG = temp[i][1] - omega[i][1] - beta * k / (2. * h * h) * (- temp[i - 1][0] - omega[i - 1][0] + 2. * temp[i][0] + 2. * omega[i][0] - temp[i + 1][0] - omega[i + 1][0]) / 2. - k * gamma / 4. * modulo * (temp[i][0] + omega[i][0]) / 2.;
    }
    if (i == N - 1){
        REAL = temp[i][0] - omega[i][0] + beta * k / (2. * h * h) * (- 2. * temp[i - 1][1] - 2. * omega[i - 1][1] + 2. * temp[i][1] + 2. * omega[i][1]) / 2. + k * gamma / 4. * modulo * (temp[i][1] + omega[i][1]) / 2.;
        IMAG = temp[i][1] - omega[i][1] - beta * k / (2. * h * h) * (- 2. * temp[i - 1][0] - 2. * omega[i - 1][0] + 2. * temp[i][0] + 2. * omega[i][0]) / 2. - k * gamma / 4. * modulo * (temp[i][0] + omega[i][0]) / 2.;
    }
    return vector <double>{REAL, IMAG};
}

vector <double> Parcial :: df1_di(int i){
    double modulo = pow((temp[i][0] + omega[i][0]), 2.) + pow((temp[i][1] + omega[i][1]), 2.);
    double dREAL = 1. + k * gamma / 4. * (temp[i][0] + omega[i][0]) * (temp[i][1] + omega[i][1]);
    double dIMAG = beta * k / (2. * h * h) + k * gamma / 4. * (temp[i][1] + omega[i][1]) * (temp[i][1] + omega[i][1]) + k * gamma * modulo / 8.;
    return vector <double>{dREAL, dIMAG};
}

vector <double> Parcial :: df2_di(int i){
    double modulo = pow((temp[i][0] + omega[i][0]), 2.) + pow((temp[i][1] + omega[i][1]), 2.);
    double dREAL = - beta * k / (2. * h * h) - k * gamma / 4. * (temp[i][0] + omega[i][0]) * (temp[i][0] + omega[i][0]) - k * gamma * modulo / 8.;
    double dIMAG = 1. - k * gamma / 2. * (temp[i][1] + omega[i][1]) * (temp[i][0] + omega[i][0]) / 2.;
    return vector <double>{dREAL, dIMAG};
}



void Parcial :: CrankNicolson( const char* filename ){

    omega.clear();

    for (int i = 0; i < N; i++){
        omega.push_back(vector <double>{f(xL + i * h).real(), f(xL + i * h).imag()}); //condiciones iniciales.
        temp.push_back(vector <double>{0, 0}); //
    }

    ofstream Solucion( filename); // El archivo solo se abre una vez
    if (!Solucion){
    cerr << "No se pudo crear el archivo \n";
    exit( 1 );
    }
    
    vector <vector <double>> b(N);

    vector <vector <vector <double>>> A(N);
    vector <vector <vector <double>>> B(2);
    vector <vector <vector <double>>> C(2); //Bloques de matrices para el jacobiano

    vector <vector <vector <double>>> L(N);
    vector <vector <vector <double>>> U(N); //Bloques de matrices de la descomposicion LU

    vector <vector <double>> z(N); //Vector que se usa para el algoritmo de Crout
    vector <vector <double>> y(N); //corrección
    double suma;

    C[0] = vector< vector <double>>{{0, - beta * k / (2 * h * h)} , {beta * k / (2 * h * h), 0}};
    C[1] = vector< vector <double>>{{0, - beta * k / (4 * h * h)} , {beta * k / (4 * h * h), 0}};

    B[0] = vector< vector <double>>{{0, - beta * k / (4 * h * h)} , {beta * k / (4 * h * h), 0}};
    B[1] = vector< vector <double>>{{0, - beta * k / (2 * h * h)} , {beta * k / (2 * h * h), 0}}; 

    bool Correcto = true;
    for (int jj = 0; jj <= m; jj++){
        temp = omega;
        Correcto = false;

        while (Correcto == false){  //Mientras la corrección de las aproximaciones no tienda a 0

            for (int i = 0; i < N; i++){
                b[i] = F(i); //Define el vector b para el sistema Jx = b que se debe resolver
            }
            for (int i = 0; i < N; i++){
                A[i] = vector< vector <double>>{df1_di(i) , df2_di(i)};
            }//Define los bloques de matrices 2x2 que forman el jacobiano tridiagonal por bloques
            
            /*
            for (int i = 0; i < N - 1; i++){
                C[i] = vector< vector <double>>{df1_dimas1(i) , df2_dimas1(i)};
            }
            for (int i = 1; i < N; i++){
                B[i] = vector< vector <double>>{df1_dimenos1(i) , df2_dimenos1(i)};
            } 
            */

           //Algoritmo de Crout por bloques de matrices 2x2
            U[0] = mat_mat_mul(inv_mat(A[0]), C[0]);
            z[0] = mat_vec_mul(inv_mat(A[0]), b[0]);

            for (int i = 1; i < N - 1; i++){
                L[i] = mat_mat_sum(A[i], mat_const_mul(mat_mat_mul(B[0], U[i - 1]), -1.));
                U[i] = mat_mat_mul(inv_mat(L[i]), C[1]);
                z[i] = mat_vec_mul(inv_mat(L[i]), vec_vec_sum(b[i], vec_const_mul(mat_vec_mul(B[0], z[i - 1]), -1.)));
            }

            L[N - 1] = mat_mat_sum(A[N- 1], mat_const_mul(mat_mat_mul(B[1], U[N - 2]), -1.));
            y[N - 1] = mat_vec_mul(inv_mat(L[N - 1]), vec_vec_sum(b[N - 1], vec_const_mul(mat_vec_mul(B[1], z[N - 2]), -1.)));


            for (int i = N - 2; i >= 0; i--){
                y[i] = vec_vec_sum(z[i], vec_const_mul(mat_vec_mul(U[i], y[i + 1]) , -1.));
            }
            suma = 0.;

            for (int i = N - 1; i >= 0; i--){
                suma = y[i][0] * y[i][0] + y[i][1] * y[i][1]; //Tamaño de la correccion que se le tuvo que hacer a la solución.
            }
            cout << "Modulo de la correccion: "<< suma << "\n";

            for (int i = N - 1; i >= 0; i--){
                temp[i][0] = temp[i][0] - y[i][0];
                temp[i][1] = temp[i][1] - y[i][1];
            }
            if (suma <= TOL){
                Correcto = true;
                cout << "Converge. Esto corresponde al paso temporal " << jj + 1 << "\n";
            }
        }
        omega = temp;

        //Solucion << jj * k;
        for (int i = 0; i < Get_N(); i++){
            if (omega[i][1] > 0){
                Solucion  << setw(9) << setprecision(8) << omega[i][0] << "+" << omega[i][1] << "i" << ","; 
            }
            else{
                Solucion << setw(9) << setprecision(8) << omega[i][0] << setw(9) << setprecision(8) << omega[i][1] << "i" << ","; 
            }
                
         }
        Solucion << "\n";
    }
    Solucion.close();

}