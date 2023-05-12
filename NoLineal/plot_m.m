filename = "./Duffin oscilator.csv";

M = table2array(readtable(filename, 'ReadVariableNames', false));

plot(M(:,1),M(:,2));
