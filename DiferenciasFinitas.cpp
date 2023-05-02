#include <iostream>
#include <cmath>
#include <iomanip>


#include "DiferenciasFinitas.h" 

DiferenciasFinitas :: DiferenciasFinitas(float x0, float xf, float y0, float yf, int N){
    Set_x0(x0);
    Set_xf(xf);
    Set_y0(y0);
    Set_yf(yf);
    Set_N(N);
}

void DiferenciasFinitas :: Set_x0(float a){
    x0 = a;
}

void DiferenciasFinitas :: Set_xf(float b){
    xf = b;
}

void DiferenciasFinitas :: Set_y0(float alpha){
    y0 = alpha;
}

void DiferenciasFinitas :: Set_yf(float beta){
    yf = beta;
}

void DiferenciasFinitas :: Set_N(int n){
    N = (n >= 2) ? n : 2;
}

float DiferenciasFinitas :: Get_x0() const{
    return x0;
}

float DiferenciasFinitas :: Get_xf() const{
    return xf;
}

float DiferenciasFinitas :: Get_y0() const{
    return y0;
}

float DiferenciasFinitas :: Get_yf() const{
    return yf;
}

int DiferenciasFinitas :: Get_N() const{
    return N;
}

int DiferenciasFinitas :: Get_MaxIter() const{
    return MaxIter;
}

std :: vector <float> DiferenciasFinitas :: Get_omega() const{
    return omega_i;
}

std :: vector <float> DiferenciasFinitas :: Get_x() const{
    return x_i;
};

void DiferenciasFinitas :: Imprimir_x() const{
    std :: cout << "Antes \n";
    std :: vector <float> x_i = Get_x();
    std :: cout << "Despues \n";
    for (float i: x_i){
        std :: cout << i << "\n";
    }
};

void DiferenciasFinitas :: Imprimir_omega() const{
    std :: cout << "Antes \n";
    std :: vector <float> omega_i = Get_omega();
    std :: cout << "Despues \n";
    for (float i: omega_i){
        std :: cout << i << "\n";
    }
};

const float DiferenciasFinitas :: p(float x) const{
    return -2 / x;
}

const float DiferenciasFinitas :: q(float x) const{
    return 2 / (x * x);
}

const float DiferenciasFinitas :: r(float x) const{
    return sin(log(x)) / (x * x);
}

const float DiferenciasFinitas :: f(float x, float y, float yp) const{
    return float(1) / 8 * (32 + 2 * x * x * x - y * yp);
}

void DiferenciasFinitas :: Set_x() {
    x_i.clear();
    for (int i = 0; i <= N + 1; i++){
        x_i.push_back(x0 + (xf - x0) / (N + 1) * i);
    }
};

void DiferenciasFinitas :: Solve_Lineal(){
    Lineal(x0, xf, y0, yf, N);
    Set_x();
}

void DiferenciasFinitas :: Solve_NoLineal(){
    No_Lineal(x0, xf, y0, yf, N);
    Set_x();
}


void DiferenciasFinitas :: Lineal(const float x0, const float xf, const float y0, const float yf, const int N){
    omega_i.clear();
    for (int i = 0; i <= N + 1; i++){
        omega_i.push_back(0);
    }
    
    float h = (xf - x0) / (N + 1); //Divide el intervalo en secciones
    float x = x0 + h;
    
    std :: vector <float> a_i(N + 2);   a_i[1] = 2 + h * h * q(x);
    std :: vector <float> b_i(N + 2);   b_i[1] = -1 + h / 2 * p(x);
    std :: vector <float> c_i(N + 2);
    std :: vector <float> d_i(N + 2);   d_i[1] = -h * h * r(x) + (1 + (h / 2) * p(x)) * y0;

    for (size_t i = 2; i <= N - 1; i++){
        x = x0 + i * h;
        a_i[i] = 2 + h * h * q(x);
        b_i[i] = -1 + (h / 2) * p(x);
        c_i[i] = -1 - (h / 2) * p(x);
        d_i[i] = -h * h * r(x);
    }

    x = xf - h;
    a_i[N] = 2 + h * h * q(x);
    c_i[N] = -1 - (h / 2) * p(x);
    d_i[N] = -h * h * r(x) + (1 - (h / 2) * p(x)) * yf;

    std :: vector <float> l_i(N + 2);   l_i[1] = a_i[1];
    std :: vector <float> u_i(N + 2);   u_i[1] = b_i[1] / a_i[1];
    std :: vector <float> z_i(N + 2);   z_i[1] = d_i[1] / l_i[1];

    for (size_t i = 2; i <= N - 1; i++){
        l_i[i] = a_i[i] - c_i[i] * u_i[i - 1];
        u_i[i] = b_i[i] / l_i[i];
        z_i[i] = (d_i[i] - c_i[i] * z_i[i - 1]) / l_i[i];
    }

    l_i[N] = a_i[N] - c_i[N] * u_i[N - 1];
    z_i[N] = (d_i[N] - c_i[N] * z_i[N - 1]) / l_i[N];

    omega_i[0] = y0;   omega_i[N + 1] = yf;    omega_i[N] = z_i[N];

    for (size_t i = N - 1; i >= 1; i--){
        omega_i[i] = z_i[i] - u_i[i] * omega_i[i + 1];
    }    
}

void DiferenciasFinitas :: No_Lineal(const float x0, const float xf, const float y0, const float yf, const int N, const float TOL, const int MaxIter){
    float h = (xf - x0) / (N + 1); //Divide el intervalo en secciones

    omega_i.clear();
    for (int i = 0; i <= N + 1; i++){
        omega_i.push_back(0);
    }
    omega_i[0] = y0;   omega_i[N + 1] = yf;

    for (size_t i = 1; i <= N; i++){
        omega_i[i] = y0 + i * (yf - y0) / (xf - x0) * h;
    }
    std :: vector <float> a_i(N + 2);
    std :: vector <float> b_i(N + 2);
    std :: vector <float> c_i(N + 2);
    std :: vector <float> d_i(N + 2);

    std :: vector <float> l_i(N + 2);
    std :: vector <float> u_i(N + 2);
    std :: vector <float> v_i(N + 2);
    std :: vector <float> z_i(N + 2);

    size_t k = 1;

    bool Completado = false;

    while ((!Completado) && (k <= MaxIter)){
        float x = x0 + h;
        float t = (omega_i[2] - y0) / (2 * h);

        a_i[1] = 2 + h * h * (f(x, omega_i[1] + h, t) - f(x, omega_i[1], t)) / h; 
        b_i[1] = -1 + (h / 2) * (f(x, omega_i[1], t + h) - f(x, omega_i[1], t)) / h; 
        d_i[1] = -(2 * omega_i[1] - omega_i[2] - y0 + h * h * f(x, omega_i[1], t));

        for (size_t i = 2; i <= N - 1; i++){
            x = x0 + i * h;
            t = (omega_i[i + 1] - omega_i[i - 1]) / (2 * h);
            a_i[i] = 2 + h * h * (f(x, omega_i[i] + h, t) - f(x, omega_i[i], t)) / h;
            b_i[i] = -1 + (h / 2) * (f(x, omega_i[i], t + h) - f(x, omega_i[i], t)) / h; 
            c_i[i] = -1 - (h / 2) * (f(x, omega_i[i], t + h) - f(x, omega_i[i], t)) / h; 
            d_i[i] = -(2 * omega_i[i] - omega_i[i + 1] - omega_i[i - 1] + h * h * f(x, omega_i[i], t));
        }

        x = xf - h;
        t = (yf - omega_i[N - 1]) / (2 * h);
        a_i[N] = 2 + h * h * (f(x, omega_i[N] + h, t) - f(x, omega_i[N], t)) / h;
        c_i[N] = -1 - (h / 2) * (f(x, omega_i[N], t + h) - f(x, omega_i[N], t)) / h;
        d_i[N] = -(2 * omega_i[N] - omega_i[N - 1] - yf + h * h * f(x, omega_i[N], t));

        l_i[1] = a_i[1];
        u_i[1] = b_i[1] / a_i[1];
        z_i[1] = d_i[1] / l_i[1];

        for (size_t i = 2; i <= N - 1; i++){
            l_i[i] = a_i[i] - c_i[i] * u_i[i - 1];
            u_i[i] = b_i[i] / l_i[i];
            z_i[i] = (d_i[i] - c_i[i] * z_i[i - 1]) / l_i[i];
        }

        l_i[N] = a_i[N] - c_i[N] * u_i[N - 1];
        z_i[N] = (d_i[N] - c_i[N] * z_i[N - 1]) / l_i[N];

        v_i[N] = z_i[N];
        omega_i[N] = omega_i[N] + v_i[N];

        for (size_t i = N - 1; i >= 1; i--){
            v_i[i] = z_i[i] - u_i[i] * v_i[i + 1];
            omega_i[i] = omega_i[i] + v_i[i];
        }

        float suma_v = 0;
        for (auto& n : v_i){
            suma_v += n * n;
        }
        
        if (suma_v <= TOL){
            bool Completado = true;
        }
        k++;
    }
    if (!Completado){
    std :: cout << "Se excedio el maximo numero de iteraciones \n";
    }
}