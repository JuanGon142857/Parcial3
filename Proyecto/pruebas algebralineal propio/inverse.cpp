#include <iostream>
#include <vector>

using namespace std;

vector <vector <double>> inv_mat(vector <vector <double>>); //Saca inverso de matriz 2x2

vector <vector <double>> mat_mat_mul(vector <vector <double>>,vector <vector <double>>);

vector <double> mat_vec_mul(vector <vector <double>>,vector <double>);

vector<vector <double>> mat_mat_sum(vector <vector <double>>, vector <vector <double>>);
vector<vector <double>> mat_const_mul(vector <vector <double>>, double);

vector <double> vec_vec_sum(vector <double>, vector <double>);
vector <double> vec_const_mul(vector <double>, double);

//vector <vector <double >> A = {{1, 0.21},{0.12 ,-5}}; 

vector <double> x ={5,6};

vector <vector <double>> b = {{1, 2},{3,4},{5,6}};
vector <vector <double>> y = b;

vector <vector <vector <double >>> A = {{{12, 25},{33 , 48}},{{41, -5},{42 , -5}}, {{5, 13},{2 , 6}}}; 
vector <vector <vector <double >>> B = {{{0, 0},{0 , 0}},{{-2, -46},{36 , -13}}, {{28, 17},{-3 , 5}}};
vector <vector <vector <double >>> C = {{{13, 1},{43 , -5}},{{29, 54},{45 , -6}}, {{0, 0},{0 , 0}}};

int main(){
    vector <vector <vector <double >>> L(3);
    vector <vector <vector <double >>> U(3);
    vector <vector <double>> z = {{1, 2},{3,4},{5,6}};
    
    L[0] = A[0];
    U[0] = mat_mat_mul(inv_mat(L[0]), C[0]);
    z[0] = mat_vec_mul(inv_mat(L[0]), b[0]);

    for (int i = 1; i < L.size() - 1; i++){
        L[i] = mat_mat_sum(A[i], mat_const_mul(mat_mat_mul(B[i], U[i - 1]), -1.));
        U[i] = mat_mat_mul(inv_mat(L[i]), C[i]);
        z[i] = mat_vec_mul(inv_mat(L[i]), vec_vec_sum(b[i], vec_const_mul(mat_vec_mul(B[i], z[i - 1]), -1.)));
    }

    L[L.size() - 1] = mat_mat_sum(A[L.size() - 1], mat_const_mul(mat_mat_mul(B[L.size() - 1], U[L.size() - 2]), -1.));
    z[L.size() - 1] = mat_vec_mul(inv_mat(L[L.size() - 1]), vec_vec_sum(b[L.size() - 1], vec_const_mul(mat_vec_mul(B[L.size() - 1], z[L.size() - 2]), -1.)));

    y[L.size() - 1] = z[L.size() - 1];

    for (int i = L.size() - 2; i >=0; i--){
        y[i] = vec_vec_sum(z[i], vec_const_mul(mat_vec_mul(U[i], y[i + 1]) , -1.));
    }
    
    cout << "L3:  \n";
    for (vector <double> i : L[2]){
        for (double j: i){
            cout << j << ", ";
        }
        cout << "\n"; 
    }

    cout << "indice : " << L.size() -1 << "\n";

    cout  << "y3 : \n";
    for (int i = 0; i <= 1; i++){
        cout << y[2][i] << "\n";
    }

    /*
    vector <double> dosx = mat_vec_mul(A, x);

    cout << "A:  \n";
    for (vector <double> i : A){
        for (double j: i){
            cout << j << ", ";
        }
        cout << "\n";
    }

    vector <vector <double>>inv = inv_mat(A);
    cout << "inv:  \n";
    for (vector <double> i : inv){
        for (double j: i){
            cout << j << ", ";
        }
        cout << "\n";
    }

    vector <vector <double>>identity = mat_mat_mul( A, inv);
    cout << "I:  \n";
    for (vector <double> i : identity){
        for (double j: i){
            cout << j << ", ";
        }
        cout << "\n";
    }

    cout  << "Ax : \n";
    for (int i = 0; i <= 1; i++){
        cout << dosx[i] << "\n";
    }*/

    return 0;
}


vector <vector <double>> inv_mat(vector <vector <double>> A){
    double det = A[0][0] * A[1][1] - A[1][0] * A[0][1];
    if (A.size() == 2 && A[0].size() == 2){
        if (det != 0){
            vector <vector <double>> inverse(2,vector<double> (2));
            inverse[0][0] = 1 / det * A[1][1];
            inverse[1][1] = 1 / det * A[0][0];
            inverse[0][1] = -1 / det * A[0][1];
            inverse[1][0] = -1 / det * A[1][0];
            return inverse;
        }
        else{
            cout << "Esta matriz no es invertible \n";
        }
    }
    else{
        cout << "Se ingresó una matriz que no es 2x2 \n";
    }
    return vector<vector <double>>{0};
    
}; //Saca inverso de matriz 2x2

vector <vector <double>> mat_mat_mul(vector <vector <double>> A,vector <vector <double>> B){
    if (A[0].size() == B.size()){
        vector <vector <double>> C(A.size(), vector<double> (B[0].size()));
        for (int i = 0; i < B[0].size(); i++){
            for (int j = 0; j < A.size(); j++){
                double suma = 0;
                for (int k = 0; k < A[0].size(); k++){
                    suma = suma + A[i][k] * B[k][j];
                } 
                C[i][j] = suma;
            }
        }
        return C;
    }
    else{
        cout << "Estas matrices no pueden multiplicarse entre s] \n";
    }
    return vector<vector <double>>{0};
}

vector <double> mat_vec_mul(vector <vector <double>> M, vector <double> v){
    if (M[0].size() == v.size()){
        vector <double> vp(M.size());
        for (int i = 0; i < v.size(); i++){
            for (int j = 0; j < M.size(); j++){
                double suma = 0;
                for (int k = 0; k < M[0].size(); k++){
                    suma = suma + M[i][k] * v[k];
                } 
                vp[i] = suma;
            }
        }
        return vp;
    }
    else{
        cout << "Estas matrices no pueden multiplicarse entre sì \n";
    }
    return vector<double>{0};
}

vector<vector <double>> mat_mat_sum(vector <vector <double>> A, vector <vector <double>> B){
    if ((A.size() == B.size()) && (A[0].size() == B[0].size())){
        vector <vector <double>> C = A;
        for (int i = 0; i < A.size(); i++){
            for (int j = 0; j < A[0].size(); j++){
                C[i][j] = A[i][j] + B[i][j];
            }
        }
        return C;
    }
    else{
        cout << "Estas matrices no pueden sumarse entre sì \n";
    }
    return vector<vector <double>>{0};
}

vector<vector <double>> mat_const_mul(vector <vector <double>> A, double c){
    vector <vector <double>> cA = A;
    for (int i = 0; i < A.size(); i++){
        for (int j = 0; j < A[0].size(); j++){
            cA[i][j] = A[i][j] * c;
        }
    }
    return cA;
}

vector <double> vec_vec_sum(vector <double> A, vector <double> B){
    if ((A.size() == B.size())){
        vector <double> C = A;
        for (int i = 0; i < A.size(); i++){
            C[i] = A[i] + B[i];
        }
        return C;
    }
    else{
        cout << "Estas matrices no pueden sumarse entre sì \n";
    }
    return vector <double>{0};
}

vector <double> vec_const_mul(vector <double> A, double c){
    vector <double> cA = A;
    for (int i = 0; i < A.size(); i++){
        cA[i] = A[i] * c;
    }
    return cA;
}