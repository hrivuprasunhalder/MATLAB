lambda_p= 1e-5;
beta= [0.000231,0.00153,0.001372,.00276,.000805,.000294];
lambda=[0.0124,0.0305,0.1115,0.301,1.138,3.01 ];
s = linspace(.001,100,1000);%reactor period

%inhour equation
rho= s.*lambda_p;
for i= 1:6
    rho=rho + s.*(beta(i)./(s+lambda(i)));
end
rho1=rho.*1e5;

figure;
plot(1./s, rho1);
set (gca,'XScale','log');
xlabel('Reactor Period T=1/s(s)');
ylabel('Reactivity (PCM)');
title('Inhour Equation Plot');
grid on;