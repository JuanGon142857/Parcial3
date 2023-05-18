#ifndef NO_LINEAL_H
#define NO_LINEAL_H

#include "../Lineal/Lineal.h"

class No_Lineal : public Lineal
{
    public:
        No_Lineal(double = 0, double = 0, double = 0, double = 0, int = 2, int = 10, double = 1E-8); // constructor predeterminado

        void Set_MaxIter( int );
        void Set_TOL ( double );
        int Get_MaxIter() const;
        double Get_TOL () const;

        const double f( double, double, double );

        void Solve();
    private:
        
        int MaxIter;
        double TOL;

};

#endif
