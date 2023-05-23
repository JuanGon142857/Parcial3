#include <vector>
#include <complex>

#ifndef PARCIAL_H
#define PARCIAL_H

class Parcial{
    public:
        Parcial(double = -30, double = 70, double = -1, double = 2./3, double = 1, int = 10, int = 10, double = 1e-4); //constructor predeterminado

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

        std :: complex <double> f(double);

        void CrankNicolson( const char* );
        void Set_X();

        std :: vector <std :: complex <double>> Get_omega() const;
        std :: vector <double> Get_X() const;

    private:

        double TOL;

        double xL; double xR;

        double l;   int N;  double h; //Parametros para la discretizacion del espacio: longitud del intervalo, numero de puntos y tamaño del paso          
        double T;   int m;  double k; //Parametros para la discretizacion del tiempo: longitud del intervalo, numero de puntos y tamaño del paso    

        double beta;
        double gamma;
        

        std :: vector <double> X;      //Puntos donde se evalua la funcion

        std :: vector <std :: vector <double>> omega;   //Aquì se guarda la parte real "v" y compleja "w" de las aproximaciones obtenidas
        //Estas partes no se puede agrupar en un vector de numeros complejos porqué parte del proceso incluye calcular el jacobiano respecto a estas variables
        //Calcular el jacobiano respecto a un número complejo resulta en un término cuya derivada no existe en los complejos


        std :: vector <std :: vector <double>> inv_mat(std :: vector <std :: vector <double>>); //Saca inverso de matriz 2x2

        std :: vector< std :: vector <double>> mat_mat_sum( std :: vector < std :: vector <double>>, std :: vector < std :: vector <double>>);
        std :: vector< std :: vector <double>> mat_const_mul( std :: vector < std :: vector <double>>, double);

        std :: vector <double> vec_vec_sum( std :: vector <double>, std :: vector <double>);
        std :: vector <double> vec_const_mul( std :: vector <double>, double);

        std :: vector <double> F(std :: vector <std :: vector <double>>, std :: vector <std :: vector <double>>, int);

        std :: vector <double> df1_di(std :: vector <std :: vector <double>>, std :: vector <std :: vector <double>>, int);
        std :: vector <double> df2_di(std :: vector <std :: vector <double>>, std :: vector <std :: vector <double>>, int);

        std :: vector <double> df1_dimenos1(std :: vector <std :: vector <double>>, std :: vector <std :: vector <double>>, int);
        std :: vector <double> df2_dimenos1(std :: vector <std :: vector <double>>, std :: vector <std :: vector <double>>, int);

        std :: vector <double> df1_dimas1(std :: vector <std :: vector <double>>, std :: vector <std :: vector <double>>, int);
        std :: vector <double> df2_dimas1(std :: vector <std :: vector <double>>, std :: vector <std :: vector <double>>, int);

        std :: vector <std :: vector <double>> mat_mat_mul(std :: vector <std :: vector <double>>, std :: vector <std :: vector <double>>);  //multiplicacion entre matrices
        std :: vector <double> mat_vec_mul(std :: vector <std :: vector <double>>, std :: vector <double>); //multiplicacion de matriz por vector

        //void Solution2csv( const char* );
};

#endif