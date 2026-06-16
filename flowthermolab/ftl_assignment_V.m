% Plot the following multi-variable function : z = xe^-(x^2 + y^2) and also plot its gradient 
% Use contour and quiver plot for this task 

clc; clear all; close all;

[X, Y] = meshgrid(-2:0.1:2, -2:0.1:2);

Z = X .* exp(-X.^2 - Y.^2);
[pX, pY] = gradient(Z, 0.1, 0.1);


figure;
hold on;
contour(X, Y, Z, 20);
quiver(X, Y, pX, pY, 'Color', 'r', 'AutoScaleFactor', 1.5);
grid on;
xlabel('X-axis');
ylabel('Y-axis');
title('Contour and Gradient(quiver) plot of z = x \cdot e^{-x^2 - y^2}');
legend('Function Contours', 'Gradient Vectors', 'Location', 'best');
axis equal;
hold off;