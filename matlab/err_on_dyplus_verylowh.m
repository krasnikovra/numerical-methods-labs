close all ; clear all; clc;
%set(0,'DefaultTextInterpreter', 'latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
%set(groot, 'defaultLegendInterpreter','latex');
%set(0,'DefaultAxesFontSize',12);
%set(0,'DefaultTextFontSize',12);

csvEuler = readmatrix("../csv/err_on_dyplus_Euler_verylowh.csv");
csvAdams = readmatrix("../csv/err_on_dyplus_Adams_verylowh.csv");

dyEuler = csvEuler(:,1);
dyAdams = csvAdams(:,1);
errEuler = csvEuler(:,2);
errAdams = csvAdams(:,2);

figure
loglog(dyEuler, errEuler, 'b',...
       dyEuler, 2.364e-08 * ones(size(dyEuler)), 'b--',...
       dyAdams, errAdams, 'r',...
       dyAdams, 4.181e-08 * ones(size(dyAdams)), 'r--')
title('Возмущение начального условия, шаг 10^{-5}', 'interpreter', 'tex')
xlabel("$$+\Delta y_0$$")
ylabel("$$\max\limits_{i=\overline{0,10^5}}|f(x_i)-y_i|$$")
legend("Мод. метод Эйлера", "Мод. метод Эйлера, без возмущения", "Метод Адамса", "Метод Адамса, без возмущения",...
       'location', 'northwest')
grid on

print -depsc ../latex/img/err_on_dyplus_verylowh.eps