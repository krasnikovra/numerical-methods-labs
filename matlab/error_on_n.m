csvData = readmatrix("../csv/error_on_n.csv");
n = csvData(1,:);
err = csvData(2,:);

figure
semilogy(n, err)
title("Error on points count")
xlabel("Points count")
ylabel("Max grid midpoints error")
grid on
