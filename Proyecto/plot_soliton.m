a = 1;
c = 1;
beta = -1;
gamma = 2./3;

x = -30 : 0.05 : 70;

phi = sqrt(2 * a / gamma) * sech(sqrt(2 * a / gamma) * x);

plot (x, abs(phi))