csvData = readmatrix("../csv/error_on_h.csv");
h = csvData(:,1);
err = csvData(:,2);

figure
loglog(h, err, ...
       h, 0.079 * h.^4, '--') % approx const is 0.079
title("Зависимость фактической ошибки от длины отрезка разбиения")
xlabel("h")
ylabel("Фактическая ошибка")
legend("Ошибка", "y=Ch^4, C=0.079", 'location', 'northwest')
grid on