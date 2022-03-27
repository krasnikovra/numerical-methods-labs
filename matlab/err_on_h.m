csvData = readmatrix("../csv/err_on_h.csv");
h = csvData(:,1);
err = csvData(:,2);

C = 228.98;

figure
loglog(h, err,...
       h, C*h.^2, '--')
title("Зависимость фактической ошибки от h")
xlabel("h")
ylabel("Фактическая ошибка")
legend("Зависимость", "y=Ch^2, C=228.98")
grid on