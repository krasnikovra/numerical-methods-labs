csvData = readmatrix("../csv/error_on_alpha.csv");
alpha = csvData(1,:);
err = csvData(2,:);

figure
semilogy(alpha, err)
title("Error on grid param for 5 points grid")
grid on