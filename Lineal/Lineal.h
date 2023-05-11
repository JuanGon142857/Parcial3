#include <vector>

#ifndef LINEAL_H
#define LINEAL_H

class Lineal{
    public:
        Lineal(float = 0, float = 0, float = 0, float = 0, int = 2); //constructor predeterminado

        void Set_x0( float );
        void Set_xf( float );
        void Set_y0( float );
        void Set_yf( float );
        void Set_N( int );

        float Get_x0() const;
        float Get_xf() const;
        float Get_y0() const;
        float Get_yf() const;
        int Get_N() const;

        const float p( float );
        const float q( float );
        const float r( float );

        void Solve();
        void Set_x();

        std :: vector <float> Get_omega() const;
        std :: vector <float> Get_x() const;

        void Solution2csv( const char* );

    protected:

        float x0;       float xf;       //limites del intervalo
        float y0;       float yf;       //funcion evaluada en los limites del intervalo 
        float N;

        std :: vector <float> x_i;      //Puntos donde se evalua la funcion
        std :: vector <float> omega_i;  //Aproximaciones obtenidas
};

#endif