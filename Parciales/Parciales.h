#include <vector>

#ifndef PARCIAL_H
#define PARCIAL_H

class Parcial{
    public:
        Parcial(float = 1, float = 1, int = 10, int = 10, float = 0); //constructor predeterminado

        void Set_l( float );
        void Set_T( float );
        void Set_m( int );
        void Set_N( int );
        void Set_alpha( float ); 

        float Get_l() const;
        float Get_T() const;
        int Get_m() const;
        int Get_N() const;
        float Get_alpha() const;

        float f(float);

        void CrankNicolson( const char*);
        void Set_x();

        std :: vector <float> Get_omega() const;
        std :: vector <float> Get_x() const;

    protected:

        float l;        float T;//limites del intervalo
        int N;          int m;
        float alpha;

        std :: vector <float> x_i;      //Puntos donde se evalua la funcion
        std :: vector <float> omega_i;  //Aproximaciones obtenidas para cada t

        void Solution2csv( const char* );
};

#endif