csvData = readmatrix("../csv/error_on_h.csv");
h = csvData(:,1);
err = csvData(:,2);

figure
loglog(h, err)
title("Зависимость фактической ошибки от длины отрезка разбиения")
xlabel("H")
ylabel("Фактическая ошибка")
legend("Ошибка", 'location', 'northwest')
grid on