csvData = readmatrix("../csv/err_on_dy.csv");
dy = csvData(:,1);
err = csvData(:,2);

figure
plot(dy, err)
title("Зависимость фактической ошибки от величины возмущения y(a)")
xlabel("\delta y(a)")
ylabel("Фактическая ошибка")
grid on