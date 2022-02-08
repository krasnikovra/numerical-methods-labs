f = @(x) 1./tan(x)-x;
csvValData = readmatrix("../csv/values.csv");
csvGridData = readmatrix("../csv/grids.csv");
n = csvGridData(1,:);
x = csvValData(1,:);
y1 = csvValData(2,:);
y2 = csvValData(3,:);
y3 = csvValData(4,:);
grid1 = csvGridData(2,:);
grid2 = csvGridData(3,:);
grid3 = csvGridData(4,:);

figure
plot(x, f(x), "--",...
     x, y1, "r",...
     x, y2, "m",...
     x, y3, "g",...
     grid1, f(grid1), "r*",...
     grid2, f(grid2), "mo",...
     grid3, f(grid3), "gs")
legend("f(x)=ctgx-x", sprintf("%d points", n(1)), sprintf("%d points", n(2)),sprintf("%d points", n(3)))
title("Hermite's polynomial interpolation")
xlabel("x")
ylabel("y")
grid on