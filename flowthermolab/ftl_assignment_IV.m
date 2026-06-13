clc;
clear all;

x = 0:0.1:7;

figure;
plot(x, sin(x), 'DisplayName', 'Sine'); hold on;
plot(x, cos(x), 'DisplayName', 'Cosine');
plot(x, tan(x), 'DisplayName', 'Tangent');
hold off; grid on;
ylim([-2, 2]);
title("sin,cos,tan behaviour");
xlabel('Angle (radians)');
ylabel('Value');
legend('Location','northeast');