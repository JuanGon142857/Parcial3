
filename = "./One Soliton.csv";

M = readtable(filename, 'ReadVariableNames', false);
M(:,1) = [];
M = table2array(M);

x = -30 : 0.05 : 70;
y = 0 :  70 * 3 / 50 / 301 : 70 * 3 / 50;

figure(1)
ax = gca;
imagesc(x, y, M)
xlabel("posición",'FontSize',25)
ylabel("tiempo",'FontSize',25)
ax.XAxis.FontSize = 25;
ax.YAxis.FontSize = 25;