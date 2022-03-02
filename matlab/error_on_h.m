csvData = readmatrix("../csv/error_on_h.csv");
h = csvData(:,1);
err = csvData(:,2);

figure
subplot(2,1,1)
loglog(h, err)
title("Зависимость фактической ошибки от длины отрезка разбиения")
xlabel("h")
ylabel("Фактическая ошибка")
grid on

csvData = readmatrix("../csv/const_approx.csv");
h = csvData(:,1);
err = csvData(:,2);

subplot(2,1,2)
semilogx(h, err)
title("Вычисление константы")
xlabel("h")
ylabel("Фактическая ошибка / h^4")
grid on