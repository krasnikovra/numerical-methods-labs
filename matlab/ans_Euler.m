close all ; clear all; clc;
set(0,'DefaultTextInterpreter', 'latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');
%set(0,'DefaultAxesFontSize',12);
%set(0,'DefaultTextFontSize',12);

f = @(x)x.^4.*(log(x)+1).^2;
x = linspace(1,2,1000);

csv1 = readmatrix("../csv/ans_Euler1.csv");
csv2 = readmatrix("../csv/ans_Euler2.csv");

x1 = csv1(:,1);
y1 = csv1(:,2);
x2 = csv2(:,1);
y2 = csv2(:,2);

figure
plot(x, f(x), ...
     x1, y1, '--*', ...
     x2, y2, '--*')
title("Решения, модифицированный метод Эйлера")
xlabel('$$x$$')
ylabel('$$y$$')
legend("Точное решение", "6 точек (h=0.2)", "11 точек (h=0.1)",...
       'location', 'northwest')
grid on

print -depsc ../latex/img/ans_Euler.eps

figure 
plot(x1, abs(f(x1)-y1), 'r--*',...
     x2, abs(f(x2)-y2), 'y--*')
title("Ошибки, модифицированный метод Эйлера")
xlabel('$$x$$')
ylabel('$$f(x_i)-y_i$$')
legend("6 точек (h=0.2)", "11 точек (h=0.1)",...
       'location', 'northwest')
grid on

print -depsc ../latex/img/err_Euler.eps