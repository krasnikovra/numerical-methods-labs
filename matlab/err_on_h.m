csvData = readmatrix("../csv/err_on_h.csv");
h = csvData(:,1);
err = csvData(:,2);

figure
loglog(h, err,...
       h, h.^2, '--')
title("Зависимость фактической ошибки от h")
xlabel("h")
ylabel("Фактическая ошибка")
legend("Зависимость", "y=h^2", 'Location', 'northwest')
grid on