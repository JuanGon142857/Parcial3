#include <vector>

#ifndef DIFERENCIASFINITAS_H
#define DIFERENCIASFINITAS_H

class DiferenciasFinitas{
    public:
        DiferenciasFinitas(float = 0, float = 0, float = 0, float = 0, int N = 2);

        void Set_x0( float );
        void Set_xf( float );
        void Set_y0( float );
        void Set_yf( float );
        void Set_N( int );
        void Set_MaxIter( int );

        const float p( float ) const;
        const float q( float ) const;
        const float r( float ) const;

        const float f(float, float, float) const;

        float Get_x0() const;
        float Get_xf() const;
        float Get_y0() const;
        float Get_yf() const;
        int Get_N() const;
        int Get_MaxIter() const;

        void Solve_Lineal();
        void Solve_NoLineal();
        void Set_x();

        void Imprimir_x() const;
        void Imprimir_omega() const;

        std :: vector <float> Get_x() const;
        std :: vector <float> Get_omega() const;

    private:

        float x0;       float xf;       //limites del intervalo
        float y0;       float yf;       //funcion evaluada en los limites del intervalo 
        float N;        float MaxIter;

        void Lineal(const float, const float, const float, const float, const int);
        void No_Lineal(const float, const float, const float, const float, const int, const float = 1E-8, const int =10); //Funciones utilitarias

        std :: vector <float> x_i;      //Puntos donde se evalua la funcion
        std :: vector <float> omega_i;  //Aproximaciones obtenidas
};

#endif