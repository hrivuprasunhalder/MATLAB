clc; clear; close all;


re = 2.8179e-13;          % cm
me = 9.109e-28;           % g
c = 2.998e10;             % cm/s
me_c2 = 0.511;            % MeV
Z = 45.798;
I = 452;              % MeV
rho = 3.67;               % g/cm^3 (for NaI)
NA = 6.022e23;            % Avogadro's number
A = 150;                  % Approx molar mass NaI
E = linspace(0,20, 500); % MeV

beta = sqrt(1 - ((me_c2./E).^2));  
Se = ((2*pi*re^2*me*c^2*NA*rho*Z/A)./(beta.^2 )).* (log((me*c^2.*(beta.^2).*E*1.602*10^-6)./(2*(I*1.602*10^-12)^2.*(1-beta.^2))) + (1-beta.^2)-((log(2)).*((beta.^2)-1+ 2.*sqrt(1-beta.^2)))+((1/8).*(1-sqrt(1-beta.^2))));
Sr = (2.*rho.*E*1.602*10^-6)./(137*me^2*c^4) .* (4.*log((2.*E.*1.602*10^-6)./(me_c2*1.602*10^-6)) - 4/3);


[~, idx]= min(Se-Sr);%idx= position in E vector where Se_vals~Sr_vals
E_equal= E(idx);
disp(['The stopping power for Beta particles in NaI is: ',num2str(E_equal)]);
                                    % note to self: print is different here

figure;hold on;
plot(E,Se,'b');
plot(E,Sr,'r');
xline(E_equal,'--k',sprintf('E= %.2f MeV\n',E_equal));
xlabel('Electron Energy(MeV)');
ylabel('Stopping power');
legend('Ionization Loss','Radiative Loss','Location','bestoutside');
title('Stopping Power vs Electron Energy(MeV)');
grid on;
hold off;
                                  