csvData = readmatrix("../csv/err_on_hpc.csv");
h = csvData(:,1);
err = csvData(:,2);

C = 81.97;

figure
loglog(h, err,...
       h, C*h.^2, '--')
title("Зависимость фактической ошибки от h")
xlabel("h")
ylabel("Фактическая ошибка")
legend("Зависимость", "y=Ch^2, C=81.97")
grid on