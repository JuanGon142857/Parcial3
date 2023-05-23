/*#include "./../src/stdafx.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "./../src/linalg.h"

using namespace alglib;


int main(int argc, char **argv)
{
    complex_2d_array a = "[[1i,-1],[1i,1]]";
    ae_int_t info;
    matinvreport rep;
    cmatrixinverse(a, info, rep);
    printf("%d\n", int(info)); // EXPECTED: 1
    printf("%s\n", a.tostring(4).c_str()); // EXPECTED: [[-0.5i,-0.5i],[-0.5,0.5]]
    printf("%.4f\n", double(rep.r1)); // EXPECTED: 0.5
    printf("%.4f\n", double(rep.rinf)); // EXPECTED: 0.5
    return 0;
}*/

/*
#include <cmath>
#include "./../src/linalg.h"
#include <complex>
#include <iomanip>
#include <iostream>
 
int main()
{
    
    using namespace std;
    std::cout << std::fixed << std::setprecision(1);
 
    std::complex<double> z1 = 1i * 1i; // imaginary unit squared
    std::cout << "i * i = " << z1 << '\n';
 
    std::complex<double> z2 = std::pow(1i, 2); // imaginary unit squared
    std::cout << "pow(i, 2) = " << z2 << '\n';
 
    const double PI = std::acos(-1); // or std::numbers::pi in C++20
    std::complex<double> z3 = std::exp(1i * PI); // Euler's formula
    std::cout << "exp(i * pi) = " << z3 << '\n';
 
    std::complex<double> z4 = 1.0 + 2i, z5 = 1.0 - 2i; // conjugates
    std::cout << "(1 + 2i) * (1 - 2i) = " << z4 * z5 << '\n';

    std::complex<double> z6 = 1.0 + 2i;
    double a = 6;
    std::cout << "Real part of -5 * (1 + 2i) = " << (1 - a) * z6.real() << '\n';

    std::complex<double> z7 = 1.0 + 3i;
    std::cout << "(1 + 3i) ** 2  = " << pow(z7, 2) << '\n';
}*/

#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <cmath>
#include <complex>
using namespace std;

const double SMALL = 1.0E-30;          // used to stop divide-by-zero
const double NEARZERO = 1.0e-10;       // helps in printing

using vec    = vector<complex <double>>;         // vector
using matrix = vector<vec>;            // matrix

// Function prototypes
void print( const string &title, const matrix &A );
matrix matmul( const matrix &A, const matrix &B );
matrix subtract( const matrix &A, const matrix &B );
matrix oppsign( matrix A );
matrix subMatrix( const matrix &A, int i1, int i2, int j1, int j2 );
matrix assembly( const matrix &A11, const matrix &A12, const matrix &A21, const matrix &A22 );
matrix inverse( const matrix &A );

//======================================================================

void print( const string &title, const matrix &A )
{
   if ( title != "" ) cout << title << '\n';

   for ( auto &row : A )
   {
      for ( auto x : row ) cout << setw( 15 ) << ( abs( x ) < NEARZERO ? 0.0 : x );
      cout << '\n';
   }
}

//======================================================================

matrix matmul( const matrix &A, const matrix &B )          // Matrix times matrix
{
   int rowsA = A.size(),   colsA = A[0].size();
   int rowsB = B.size(),   colsB = B[0].size();
   assert( colsA == rowsB );

   matrix C( rowsA, vec( colsB, 0.0 ) );
   for ( int i = 0; i < rowsA; i++ )
   {
      for ( int j = 0; j < colsB; j++ )
      {
         for ( int k = 0; k < colsA; k++ ) C[i][j] += A[i][k] * B[k][j];
      }
   }
   return C;
}

//======================================================================

matrix subtract( const matrix &A, const matrix &B )        // Subtract matrices
{
   int rows = A.size(),   cols = A[0].size();
   assert( rows == B.size() && cols == B[0].size() );

   matrix result( rows, vec( cols ) );
   for ( int i = 0; i < rows; i++ )
   {
      for ( int j = 0; j < cols; j++ ) result[i][j] = A[i][j] - B[i][j];
   }
   return result;
}

//======================================================================

matrix oppsign( matrix A )                                  // Minus matrix
{
   for ( auto &row : A )
   {
      for ( auto &e : row ) e = -e;
   }
   return A;
}

//======================================================================

matrix subMatrix( const matrix &A, int i1, int i2, int j1, int j2 )
{
   int rows = i2 - i1 + 1, cols = j2 - j1 + 1;
   matrix result( rows, vec( cols ) );
   for ( int i = i1, r = 0; i <= i2; i++, r++ )
   {
      auto it1 = A[i].begin() + j1, it2 = A[i].begin() + j2 + 1;
      copy( it1, it2, result[r].begin() );
   }
   return result;
}

//======================================================================

matrix assembly( const matrix &A11, const matrix &A12, const matrix &A21, const matrix &A22 )
{
   int k = A11.size();           
   int n = k + A22.size();
   matrix result( n, vec( n ) );

   for ( int i = 0; i < k; i++ )
   {
      copy( A11[i].begin(), A11[i].end(), result[i].begin()     );
      copy( A12[i].begin(), A12[i].end(), result[i].begin() + k );
   }

   for ( int i = k; i < n; i++ )
   {
      copy( A21[i-k].begin(), A21[i-k].end(), result[i].begin()     );
      copy( A22[i-k].begin(), A22[i-k].end(), result[i].begin() + k );
   }

   return result;
}

//======================================================================

matrix inverse( const matrix &A )
{
   int n = A.size();
   if ( n == 1 ) 
   {
      complex <double> value = A[0][0];
      if ( abs( value ) < SMALL )
      {
         cerr << "Non-invertible. Giving up.\n";
         exit( 0 );
      }
      return matrix( 1, vec( 1, 1.0 / value ) );
   }

   // Partition into four
   int k = n / 2;
   matrix A11 = subMatrix( A, 0, k - 1, 0, k - 1 );
   matrix A12 = subMatrix( A, 0, k - 1, k, n - 1 );
   matrix A21 = subMatrix( A, k, n - 1, 0, k - 1 );
   matrix A22 = subMatrix( A, k, n - 1, k, n - 1 );

   // Strassen steps
   matrix R1  = inverse( A11 );
   matrix R2  = matmul( A21, R1 );
   matrix R3  = matmul( R1, A12 );
   matrix R4  = matmul( A21, R3 );
   matrix R5  = subtract( R4, A22 );
   matrix R6  = inverse( R5 );
   matrix X12 = matmul( R3, R6 );
   matrix X21 = matmul( R6, R2 );
   matrix R7  = matmul( R3, X21 );
   matrix X11 = subtract( R1, R7 );
   matrix X22 = oppsign( R6 );

   return assembly( X11, X12, X21, X22);
}

//====================================================================== 