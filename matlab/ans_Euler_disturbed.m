set(0,'DefaultTextInterpreter', 'latex');
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');
%set(0,'DefaultAxesFontSize',12);
%set(0,'DefaultTextFontSize',12);

f = @(x)x.^4.*(log(x)+1).^2;
x = linspace(1,2,1000);

csv2 = readmatrix("../csv/ans_Euler2.csv"); % h = 0.1 no distubance
csvDist = readmatrix("../csv/ans_Euler_disturbed1.csv"); % h = 0.1 disturbed by +0.1

x2 = csv2(:,1);
y2 = csv2(:,2);
xDist = csvDist(:,1);
yDist = csvDist(:,2);

figure
plot(x, f(x), ...
     x2, y2, 'm--*', ...
     xDist, yDist, 'k--*')
title('Решение возмущенной задачи, мод. метод Эйлера, h = 0.1')
xlabel('$$x$$')
ylabel('$$y$$')
legend("Точное решение", "$$\Delta y_0 = 0$$", "$$\Delta y_0 = 0.1$$",...
       'location', 'northwest')
grid on

print -depsc ../latex/img/ans_Euler_disturbed.eps

figure
plot(x2, f(x2)-y2, 'm--*', ...
     xDist, f(xDist)-yDist, 'k--*')
title('Погрешность решения возмущенной задачи, мод. метод Эйлера, h = 0.1')
xlabel('$$x$$')
ylabel('$$f(x_i)-y_i$$')
legend("$$\Delta y_0 = 0$$", "$$\Delta y_0 = 0.1$$",...
       'location', 'northwest')
grid on

print -depsc ../latex/img/err_Euler_disturbed.eps
