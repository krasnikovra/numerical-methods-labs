csvData = readmatrix("../csv/q_on_eps.csv");
eps = csvData(:,1);
q = csvData(:,2);

figure
semilogx(eps, q)
title("Зависимость числа итераций от задаваемой точности")
xlabel("\epsilon")
ylabel("q = log_2m")
grid on