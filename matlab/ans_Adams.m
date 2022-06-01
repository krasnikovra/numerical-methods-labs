close all ; clear all; clc;
set(0,'DefaultTextInterpreter', 'latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');
%set(0,'DefaultAxesFontSize',12);
%set(0,'DefaultTextFontSize',12);

f = @(x)x.^4.*(log(x)+1).^2;
x = linspace(1,2,1000);

csv1 = readmatrix("../csv/ans_Adams1.csv");
csv2 = readmatrix("../csv/ans_Adams2.csv");

x1 = csv1(:,1);
y1 = csv1(:,2);
x2 = csv2(:,1);
y2 = csv2(:,2);

figure
plot(x, f(x), ...
     x1, y1, 'r-.s', ...
     x2, y2, 'g-.s')
title("Решения, явный метод Адамса 2-го порядка")
xlabel('$$x$$')
ylabel('$$y$$')
legend("Точное решение", "6 точек (h=0.2)", "11 точек (h=0.1)",...
       'location', 'northwest')
grid on

print -depsc ../latex/img/ans_Adams.eps