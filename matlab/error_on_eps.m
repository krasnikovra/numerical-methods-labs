csvData = readmatrix("../csv/error_on_eps.csv");
eps = csvData(:,1);
err = csvData(:,2);
err_rich = csvData(:,3);

figure
loglog(eps, err,...,
       eps, err_rich,...,
       eps, eps, "--")
title("Зависимость фактической ошибки от задаваемой точности")
xlabel("Задаваемая точность \epsilon")
ylabel("Фактическая ошибка")
legend("Ошибка", "Ошибка (с поправкой Ричардсона)", "y=x", 'location', 'northwest')
grid on