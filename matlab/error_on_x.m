f = @(x) 1./tan(x) - x;
csvValData = readmatrix("../csv/values.csv");
csvGridData = readmatrix("../csv/grids.csv");
csvErrorTheoreticData = readmatrix("../csv/error_theoretic.csv");
err_theor = csvErrorTheoreticData(1,:)
n = csvGridData(1,:); % look it the file
x = csvValData(1,:);
y1 = csvValData(2,:);
y2 = csvValData(3,:);
y3 = csvValData(4,:);
grid1 = csvGridData(2,:);
grid2 = csvGridData(3,:);
grid3 = csvGridData(4,:);

fx = f(x);

figure
plot(x, abs(y1-fx), "r",...
     x, abs(y2-fx), "m",...
     x, abs(y3-fx), "g")
xlabel("x")
ylabel("|P_n(x)-f(x)|")
title(sprintf("Fact error on x for n=%d,%d,%d polynomials", n(1),n(2),n(3)))
legend(sprintf("%d points", n(1)), sprintf("%d points", n(2)),sprintf("%d points", n(3)), "theoretical for n=3")
grid on
