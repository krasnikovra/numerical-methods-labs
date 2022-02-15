f = @(x) 1./tan(x)-x;
csvValData = readmatrix("../csv/values.csv");
n = [4 5 6];
x = csvValData(1,:);
y1 = csvValData(2,:);
y2 = csvValData(3,:);
y3 = csvValData(4,:);

fx = f(x);

figure
plot(x, abs(y1-fx), "r",...
     x, abs(y2-fx), "m",...
     x, abs(y3-fx), "g")
xlabel("x")
ylabel("|P_n(x)-f(x)|")
title(sprintf("Fact error on x for n=%d,%d,%d polynomials", n(1),n(2),n(3)))
legend(sprintf("%d points", n(1)), sprintf("%d points", n(2)),sprintf("%d points", n(3)))
grid on
