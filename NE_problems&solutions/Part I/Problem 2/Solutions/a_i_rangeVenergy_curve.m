clc;clear;close all;

E_low= linspace(0.01, 2.5, 200);
E_high= linspace(2.5, 10,200);

R_low= 412.*E_low.^(1.265-0.0954.*log(E_low));
R_high= (530.*E_high)-106;

E= [E_low E_high];
R= [R_low R_high];

% Densities of materials(g/cm^3)
% air= .001225, polyethylene= 0.94, aluminum= 2.70, iron= 7.87, lead= 11.34

% Converting range from mg/cm^2 to cm for each material
% R_mass= density*R_linear, R_linear= R_mass/density
% R_mass is in mg/cm^2, we turn it into g/cm^2 so that we can divide by the density
% So R_linear= R_mass*e-3/density

R_air= (R.*10^-3)./0.001225;
R_poly= (R.*10^-3)./0.94;
R_Al= (R.*10^-3)./2.70;
R_Fe= (R.*10^-3)./7.87;
R_Pb= (R.*10^-3)./11.34;

figure;
hold on;

plot(E, R_air, 'b','LineWidth',1.5);
plot(E, R_poly, 'g','LineWidth',1.5);
plot(E, R_Al, 'r','LineWidth',1.5);
plot(E, R_Fe, 'm','LineWidth',1.5);
plot(E, R_Pb, 'w','LineWidth',1.5);
set(gca, 'YScale','log'); % gca= get cueent axes,This is a must if curvature in the plot is required
hold off;

xlabel('Beta Particle Energy(MeV)');
ylabel('Range(cm)');
title("Range-Energy Curve for Beta Particles in Different Materials");
legend('Air','Polyethylene','Aluminum','Iron','Lead','Location','bestoutside');
grid on;
