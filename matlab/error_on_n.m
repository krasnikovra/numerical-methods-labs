csvData = readmatrix("../csv/error_on_n.csv");
n = csvData(:,1);
err = csvData(:,2);

figure
semilogy(n, err)
title("Зависимость фактической ошибки от числа отрезков разбиения в виде 2^n")
xlabel("n = log_2 числа отрезков разбиения")
ylabel("Фактическая ошибка")
grid on