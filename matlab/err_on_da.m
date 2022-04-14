csvData = readmatrix("../csv/err_on_da.csv");
da = csvData(:,1);
err = csvData(:,2);

figure
subplot(2,1,1)
loglog(da, err)
title("Зависимость фактической ошибки от величины возмущения A")
xlabel("\delta A")
ylabel("Фактическая ошибка")
grid on

csvData = readmatrix("../csv/err_on_db.csv");
db = csvData(:,1);
err = csvData(:,2);

subplot(2,1,2)
loglog(db, err)
title("Зависимость фактической ошибки от величины возмущения B")
xlabel("\delta B")
ylabel("Фактическая ошибка")
grid on