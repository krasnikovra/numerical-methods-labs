close all ; clear all; clc;
set(0,'DefaultTextInterpreter', 'latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');
%set(0,'DefaultAxesFontSize',12);
%set(0,'DefaultTextFontSize',12);

f = @(x)x.^4.*(log(x)+1).^2;

csvEuler1 = readmatrix("../csv/ans_Euler1.csv");
csvEuler2 = readmatrix("../csv/ans_Euler2.csv");
csvAdams1 = readmatrix("../csv/ans_Adams1.csv");
csvAdams2 = readmatrix("../csv/ans_Adams2.csv");

xE1 = csvEuler1(:,1);
yE1 = csvEuler1(:,2);
xE2 = csvEuler2(:,1);
yE2 = csvEuler2(:,2);

xA1 = csvAdams1(:,1);
yA1 = csvAdams1(:,2);
xA2 = csvAdams2(:,1);
yA2 = csvAdams2(:,2);

figure
plot(xE1, f(xE1)-yE1, 'm--*',...
     xE2, f(xE2)-yE2, 'k--*',...
     xA1, f(xA1)-yA1, 'r-.s',...
     xA2, f(xA2)-yA2, 'g-.s')
title('Погрешности полученных решений')
xlabel('$$x$$')
ylabel('$$f(x_i)-y_i$$')
legend('Мод. метод Эйлера, h = 0.2',...
       'Мод. метод Эйлера, h = 0.1',...
       'Метод Адамса, h = 0.2',...
       'Метод Адамса, h = 0.1',...
       'location', 'northwest')
grid on

print -depsc ../latex/img/err_on_x.eps