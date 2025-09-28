% BE/A SEMF vs EXPERIMENTAL
% Excluded eta term in SEMF due to its insignificance
clear; close all; clc;
av= 16;
as= 18;
ac= 0.72;
aa= 23.5;
ap= 11;                                                         % all in MeV
data= readtable('b_binding_energy_per_nucleon.csv');
A= data.A;
Z= data.Z;
BEAEXP= (data.B_A_EXP_MeV_ );                                   %BE/A Experimental
N= A- Z;
vol= av .*A;
surf= -as.* (A.^ (2/ 3));
coulomb= -ac.* (Z.* (Z- 1))./ (A.^ (1/ 3));
asym= -aa.* ((N- Z).^ 2)./ A;
delta= (((-1).^ Z)+ ((-1).^ N)).* (ap./ sqrt(A));                                 
BEv= vol;
BEvs= vol+ surf;
BEvsc= vol+ surf+ coulomb;
BEvsca= vol+ surf+ coulomb+ asym;
BEvscad= vol+ surf+ coulomb+ asym+ delta;

figure('Units','pixels','Position',[100 100 800 600]);hold on;box on;

plot(A,BEv./ A , '-','LineWidth',1,'DisplayName','BE/A(v)');
plot(A,BEvs./ A, '-','LineWidth',1,'DisplayName','BE/A(vs)');
plot(A,BEvsc./ A, '-','LineWidth',1,'DisplayName','BE/A(vsc)');
plot(A,BEvsca./ A, '-','LineWidth',1,'DisplayName','BE/A(vsca)');  %BE convertion to BE/A
plot(A,BEvscad./A, '-','LineWidth',2,'DisplayName','BE/A(vscad)');%The complete equation
scatter(A,BEAEXP,40,'rx','LineWidth',2,'DisplayName','B/A(exp)');

xlabel('Mass Number A','FontSize',12);
ylabel('Binding Energy per Nucleon(MeV)','FontSize',12);
title('SEMF vs Experimental data','FontSize',13);
legend('Location','bestoutside');
grid on;
hold off;




