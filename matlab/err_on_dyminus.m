close all ; clear all; clc;
set(0,'DefaultTextInterpreter', 'latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');
%set(0,'DefaultAxesFontSize',12);
%set(0,'DefaultTextFontSize',12);

csvEuler = readmatrix("../csv/err_on_dyminus_Euler.csv");
csvAdams = readmatrix("../csv/err_on_dyminus_Adams.csv");

dyEuler = csvEuler(:,1);
dyAdams = csvAdams(:,1);
errEuler = csvEuler(:,2);
errAdams = csvAdams(:,2);

figure
loglog(dyEuler, errEuler,...
       dyAdams, errAdams)
title("Возмущение начального условия")
xlabel("$$-\Delta y_0$$")
ylabel("$$\max\limits_{i=\overline{0,10}}|f(x_i)-y_i|$$")
legend("Мод. метод Эйлера", "Метод Адамса", 'location', 'northwest')
grid on

print -depsc ../latex/img/err_on_dyminus.eps