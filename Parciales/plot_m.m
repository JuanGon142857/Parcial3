filename = "./Heat equation.csv";

M = readtable(filename, 'ReadVariableNames', false);
M(:,1) = [];
M = table2array(M);

x = 1/100 : 1/100 : 1;
y = 1/100 : 1/100 : 1;

figure(1)
surf(x, y, M)
xlabel("longitud")
ylabel("tiempo")
% surf(x, y, 1 / pi ^ 2 * exp(-meshgrid(y)') .* sin(pi * meshgrid(x)))

figure(2)
plot(sum(abs(M - 1 / pi ^ 2 * exp(-meshgrid(y)') .* sin(pi * meshgrid(x)))))