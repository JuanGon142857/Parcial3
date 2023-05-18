filename = "./Heat equation.csv";

M = readtable(filename, 'ReadVariableNames', false);
M(:,1) = [];
M = table2array(M);

x = 1/100 : 1/100 : 1;
y = 1/100 : 1/100 : 1;

xinverse = 1: -1/100 : 1/100;

figure(1)
ax = gca; 
ax.FontSize = 25; 
surf(xinverse, y, M)
xlabel("posición relativa",'FontSize',25)
% xticks('FontSize',25)
ylabel("tiempo (s)",'FontSize',25)
zlabel("Diferencia de temperatura respecto al reservorio (u.a)",'FontSize',20)
ax = gca(1); 
ax.XAxis.FontSize = 25;
ax.YAxis.FontSize = 25;
ax.ZAxis.FontSize = 25;
% surf(x, y, 1 / pi ^ 2 * exp(-meshgrid(y)') .* sin(pi * meshgrid(x)))

% figure(2)
% plot(sum(abs(M - 1 / pi ^ 2 * exp(-meshgrid(y)') .* sin(pi * meshgrid(x)))))