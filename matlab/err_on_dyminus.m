close all ; clear all; clc;
%set(0,'DefaultTextInterpreter', 'latex');
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
loglog(dyEuler, errEuler, 'b',...
       dyEuler, 1.763 * ones(size(dyEuler)), 'b--',...
       dyAdams, errAdams, 'r',...
       dyAdams, 3.055 * ones(size(dyAdams)), 'r--')
title("Возмущение начального условия, шаг 0.1")
xlabel("$$-\Delta y_0$$")
ylabel("$$\max\limits_{i=\overline{0,10}}|f(x_i)-y_i|$$")
legend("Мод. метод Эйлера", "Мод. метод Эйлера, без возмущения", "Метод Адамса", "Метод Адамса, без возмущения",...
       'location', 'northwest')
grid on

print -depsc ../latex/img/err_on_dyminus.eps