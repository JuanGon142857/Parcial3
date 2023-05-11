#ifndef NO_LINEAL_H
#define NO_LINEAL_H

#include "../Lineal/Lineal.h"

class No_Lineal : public Lineal
{
    public:
        No_Lineal(float = 0, float = 0, float = 0, float = 0, int = 2, int = 10, float = 1E-8); // constructor predeterminado

        void Set_MaxIter( int );
        void Set_TOL ( float );
        int Get_MaxIter() const;
        float Get_TOL () const;

        const float f( float, float, float );

        void Solve();
    private:
        
        int MaxIter;
        float TOL;

};

#endif
