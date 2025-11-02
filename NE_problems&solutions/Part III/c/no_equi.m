clear; clc; close all

% Secular equilibrium: 138Xe -> 138Cs
T1_parent = 14.08*60;
T1_daughter = 33.41*60;

lambda1 = log(2)/T1_parent;
lambda2 = log(2)/T1_daughter;

N10 = 1;
N20 = 0;

t_end = T1_daughter*10;
dt = 50;

t = 0:dt:t_end;
N1 = zeros(size(t));
N2 = zeros(size(t));
N1(1) = N10;
N2(1) = N20;

f1 = @(N1) -lambda1*N1;
f2 = @(N1, N2) lambda1*N1 - lambda2*N2;

for i = 1:length(t)-1

    %---- Parent RK4 ----
    k1 = f1(N1(i));
    k2 = f1(N1(i) + k1*dt/2);
    k3 = f1(N1(i) + k2*dt/2);
    k4 = f1(N1(i) + k3*dt);
    N1(i+1) = N1(i) + (dt/6)*(k1 + 2*k2 + 2*k3 + k4);

    %---- Daughter RK4 ----
    m1 = f2(N1(i), N2(i));
    m2 = f2(N1(i) + k1*dt/2, N2(i) + m1*dt/2);
    m3 = f2(N1(i) + k2*dt/2, N2(i) + m2*dt/2);
    m4 = f2(N1(i) + k3*dt,   N2(i) + m3*dt);
    N2(i+1) = N2(i) + (dt/6)*(m1 + 2*m2 + 2*m3 + m4);

end

% Activities

A1 = (lambda1.*N1)./2.818*10^08;
A2 = (lambda2.*N2)./2.818*10^08;


figure; hold on;
plot(t/3600, A1,'r','LineWidth',1.5,'DisplayName','138-Xe');
plot(t/3600, A2,'b','LineWidth',1.5,'DisplayName','138-Cs');
xlabel('Time (hours)');
ylabel('Activity');
title('No Equilibrium: ^{138}Xe -> ^{138}Cs');
legend('Location','bestoutside');
grid on;
hold off;
