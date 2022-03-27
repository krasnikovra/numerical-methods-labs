csvData = readmatrix("../csv/err_on_eps.csv");
eps = csvData(:,1);
err = csvData(:,2);
err1 = csvData(:,3);

figure
loglog(eps, err,...
       eps, err1,...
       eps, eps, "--")
title("Зависимость фактической ошибки от задаваемой точности")
xlabel("Задаваемая точность \epsilon")
ylabel("Фактическая ошибка")
legend("Макс. ошибка", "Ошибка во второй точке", "y=x", 'location', 'northwest')
grid on