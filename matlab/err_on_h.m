close all ; clear all; clc;
set(0,'DefaultTextInterpreter', 'latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');
%set(0,'DefaultAxesFontSize',12);
%set(0,'DefaultTextFontSize',12);

csvEuler = readmatrix("../csv/err_on_h_Euler.csv");
csvAdams = readmatrix("../csv/err_on_h_Adams.csv");
hEuler = csvEuler(:,1);
hAdams = csvAdams(:,1);
errEuler = csvEuler(:,2);
errAdams = csvAdams(:,2);

figure
loglog(hEuler, errEuler,...
       hAdams, errAdams,...
       hAdams, 419*hAdams.^2, '--',...
       hEuler, 229*hAdams.^2, '--')
title("Зависимость нормы погрешности решения от шага")
xlabel("$$h$$")
ylabel("$$\max\limits_{i=\overline{0,10}}|f(x_i)-y_i|$$")
legend("Мод. метод Эйлера", "Метод Адамса", "$y=419h^2$", "$y=229h^2$", 'location', 'Northwest')
grid on

print -depsc ../latex/img/err_on_h.eps