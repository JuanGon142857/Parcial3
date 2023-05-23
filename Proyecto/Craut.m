A = [[12, 25, 13, 1, 0 ,0]; 
    [33, 48, 43, -5, 0, 0]; 
    [-2, -46, 41, -5, 29, 54]; 
    [36, -13, 42, -5, 45, -6]; 
    [0, 0, 28, 17, 5, 13]; 
    [0, 0, -3, 5, 2, 6]];

b = [1,2,3,4,5,6]';

A1 = [[12, 25]; [33, 48]];
A2 = [[41, -5]; [42, -5]];
A3 = [[5, 13]; [2, 6]];

C1 = [[13, 1]; [43, -5]];
C2 = [[29, 54]; [45, -6]];

B2 = [[-2, -46]; [36, -13]];
B3 = [[28, 17]; [-3, 5]];

b1 = [1, 2]';
b2 = [3, 4]';
b3 = [5, 6]';

L1 = A1;
U1 = inv(L1) * C1;
z1 = inv(L1) * b1;

L2 = A2 - B2 * U1;
U2 = inv(L2) * C2;
z2 = inv(L2) * (b2 - B2 * z1);

L3 = A3 - B3 * U2;
z3 = inv(L3) * (b3 - B3 * z2);

x3 = z3;
x2 = z2 - U2 * x3;
x1 = z1 - U1 * x2;

