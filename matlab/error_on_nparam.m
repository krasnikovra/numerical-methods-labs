csvData = readmatrix("../csv/error_on_n.csv");
n = csvData(1,:);
err = csvData(2,:);

csvData = readmatrix("../csv/error_on_nparam.csv");
nparam = csvData(1,:);
errparam = csvData(2,:);

figure
semilogy(n, err, "--",...
         nparam, errparam)
title("Error on points count")
xlabel("Points count")
ylabel("Max grid midpoints error")
legend("\alpha=1", "\alpha=2")
grid on