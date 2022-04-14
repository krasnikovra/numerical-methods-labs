csvData = readmatrix("../csv/err_on_db.csv");
db = csvData(:,1);
err = csvData(:,2);

figure
loglog(db, err)
title("Зависимость фактической ошибки от величины возмущения B")
xlabel("\delta B")
ylabel("Фактическая ошибка")
grid on