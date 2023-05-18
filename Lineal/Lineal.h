#include <vector>

#ifndef LINEAL_H
#define LINEAL_H

class Lineal{
    public:
        Lineal(double = 0, double = 0, double = 0, double = 0, int = 2); //constructor predeterminado

        void Set_x0( double );
        void Set_xf( double );
        void Set_y0( double );
        void Set_yf( double );
        void Set_N( int );

        double Get_x0() const;
        double Get_xf() const;
        double Get_y0() const;
        double Get_yf() const;
        int Get_N() const;

        const double p( double );
        const double q( double );
        const double r( double );

        void Solve();
        void Set_x();

        std :: vector <double> Get_omega() const;
        std :: vector <double> Get_x() const;

        void Solution2csv( const char* );

    protected:

        double x0;       double xf;       //limites del intervalo
        double y0;       double yf;       //funcion evaluada en los limites del intervalo 
        int N;

        std :: vector <double> x_i;      //Puntos donde se evalua la funcion
        std :: vector <double> omega_i;  //Aproximaciones obtenidas
};

#endif