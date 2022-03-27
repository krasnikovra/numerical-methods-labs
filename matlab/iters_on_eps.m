csvData = readmatrix("../csv/iters_on_eps.csv");
eps = csvData(:,1);
iters = csvData(:,2);

figure
semilogx(eps, iters)
title("Зависимость числа итераций от задаваемой точности")
xlabel("\epsilon")
ylabel("Максимальное число итераций")
grid on