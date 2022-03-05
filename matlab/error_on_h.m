csvData = readmatrix("../csv/error_on_h.csv");
h = csvData(:,1);
err = csvData(:,2);
min = csvData(:,3);
max = csvData(:,4);

figure
loglog(h, err, ...
       h, 0.079 * h.^4, '--',... % approx const is 0.079
       h, max) 
title("Зависимость фактической ошибки от длины отрезка разбиения")
xlabel("H")
ylabel("Фактическая ошибка")
legend("Ошибка", "y=Ch^4, C=0.079", "Максимальная ошибка", 'location', 'northwest')
grid on