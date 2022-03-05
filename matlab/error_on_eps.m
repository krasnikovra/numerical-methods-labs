csvData = readmatrix("../csv/error_on_eps.csv");
eps = csvData(:,1);
err = csvData(:,2);
min = csvData(:,3);
max = csvData(:,4);

figure
loglog(eps, err,...
       eps, max,...
       eps, eps, "--")
title("Зависимость фактической ошибки от задаваемой точности")
xlabel("Задаваемая точность \epsilon")
ylabel("Фактическая ошибка")
legend("Ошибка", "Максимальная ошибка", "y=x", 'location', 'northwest')
grid on