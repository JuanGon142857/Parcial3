#include <vector>
#include <complex>

#ifndef PARCIAL_H
#define PARCIAL_H

class Parcial{
    public:
        Parcial(double, double, double, double, double, int, int, double, std :: complex <double> (*)(double)); //constructor predeterminado

        void Set_xL( double );
        void Set_xR( double );
        void Set_beta( double);
        void Set_gamma( double);

        void Set_T( double );
        void Set_m( int );
        void Set_N( int );
        void Set_TOL (double);

        double Get_T() const;
        int Get_m() const;
        int Get_N() const;

        void CrankNicolson( const char* );
        void Set_X();

        std :: vector <double> Get_X() const;

    private:

        double TOL;

        double xL; double xR;

        double l;   int N;  double h; //Parametros para la discretizacion del espacio: longitud del intervalo, numero de puntos y tamaño del paso          
        double T;   int m;  double k; //Parametros para la discretizacion del tiempo: longitud del intervalo, numero de puntos y tamaño del paso    

        double beta;
        double gamma;
        

        std :: vector <double> X;      //Puntos donde se evalua la funcion

        std :: vector <std :: vector <double>> omega{};   //Aquì se guarda la parte real "v" y compleja "w" de las aproximaciones obtenidas
        std :: vector <std :: vector <double>> temp{};    //Aqui se guardan las aproximaciones de cada paso mientras la solucion converge al valor buscado
        //Estas partes no se puede agrupar en un vector de numeros complejos porqué parte del proceso incluye calcular el jacobiano respecto a estas variables
        //Calcular el jacobiano respecto a un número complejo resulta en un término cuya derivada no existe en los complejos


        std :: vector <std :: vector <double>> inv_mat(std :: vector <std :: vector <double>>); //Saca inverso de matriz 2x2

        std :: vector< std :: vector <double>> mat_mat_sum( std :: vector < std :: vector <double>>, std :: vector < std :: vector <double>>);
        std :: vector< std :: vector <double>> mat_const_mul( std :: vector < std :: vector <double>>, double);

        std :: vector <double> mat_mat_sum( std :: vector <double>, std :: vector <double>);
        std :: vector <double> mat_const_mul( std :: vector <double>, double);

        std :: vector <double> F(int);

        std :: vector <double> df1_di(int);
        std :: vector <double> df2_di(int);

        std :: vector <std :: vector <double>> mat_mat_mul(std :: vector <std :: vector <double>>, std :: vector <std :: vector <double>>);  //multiplicacion entre matrices
        std :: vector <double> mat_mat_mul(std :: vector <std :: vector <double>>, std :: vector <double>); //multiplicacion de matriz por vector

        //void Solution2csv( const char* );
};

#endif