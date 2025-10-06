clc; clear; close all;
E= linspace(.01,10,400);

%atomic number_effective & approx
z_air= 7.3;
z_poly= 6;
z_Al= 13;
z_Fe= 26;
z_Pb= 82;

%Bremsstrahlung fraction= (E*Z)/700;
f_air= (E.*z_air)./700;
f_poly= (E.*z_poly)./700;
f_Al= (E.*z_Al)./700;
f_Fe= (E.*z_Fe)./700;
f_Pb= (E.*z_Pb)./700;

figure;
hold on;

plot(E, f_air, 'b','LineWidth',1.5);
plot(E, f_poly, 'g','LineWidth',1.5);
plot(E, f_Al, 'r','LineWidth',1.5);
plot(E, f_Fe, 'm','LineWidth',1.5);
plot(E, f_Pb, 'w','LineWidth',1.5);
set(gca, 'YScale','log');

hold off;

xlabel('Energy(MeV)');
ylabel('Bremsstrahlung Fraction');
title('Bremsstrahlung Fraction-Energy Curve for Beta Particles');
legend('Air','Polyethylene','Aluminum','Iron','Lead','Location','bestoutside');
grid on;
