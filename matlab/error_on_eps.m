csvData = readmatrix("../csv/error_on_eps.csv");
eps = csvData(:,1);
err = csvData(:,2);

figure
loglog(eps, err,...
       eps, eps, "--")
title("Зависимость фактической ошибки от задаваемой точности")
xlabel("Задаваемая точность \epsilon")
ylabel("Фактическая ошибка")
grid on