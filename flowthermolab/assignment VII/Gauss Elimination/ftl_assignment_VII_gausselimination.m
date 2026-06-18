% Gauss Elimination Method
clc; clear; close all;

% Co-eff matrix A
A = [2 2 -1 1;
     4 3 -1 2;
     8 5 -3 4;
     3 3 -2 2] ;

% constants b
b = [4; 6; 12; 6 ];

% Augmented matrix [A|b]
Ab = [A b];

%number of equations
n = size(A, 1);

for k = 1:n-1
    %pivoting
    [~, idx] = max(abs(Ab(k:n,k)));
    idx = idx + k - 1;
    if idx ~= k
        Ab([k idx],:) = Ab([idx k],:);
    end
    %elimination
    for i = k+1:n
        factor = Ab(i,k) / Ab(k,k);
        Ab(i,:) = Ab(i, :) - factor * Ab(k, :);
    end
end

% Back substitution
x_gauss = zeros(n,1);
for i = n:-1:1
    x_gauss(i) = (Ab(i,end) - Ab(i,1:n)*x_gauss) / Ab(i,i);
end

disp('Solution using Gaussian Elimination:');
disp(x_gauss)

%inverse method
x_inv = inv(A) * b;

disp('Solution using Matrix Inversion:');
disp(x_inv);
