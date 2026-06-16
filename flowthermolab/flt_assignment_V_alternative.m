% Plot the following multi-variable function : z = xe^-(x^2 + y^2) and also plot its gradient 
% Use contour and quiver plot for this task 

clc; clear all; close all;

[X, Y] = meshgrid(-2:0.1:2, -2:0.1:2);

Z = X .* exp(-X.^2 - Y.^2);
[pX, pY] = gradient(Z, 0.1, 0.1);

figure(1)
surf(X, Y, Z);
xlabel('X');
ylabel('Y');
zlabel('Z');

figure(2)
contour(X, Y, Z); hold on;
quiver( X, Y, pX, pY);
xlabel('X');
ylabel('Y');
hold off;
