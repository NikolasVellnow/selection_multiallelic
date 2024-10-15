% Script to play around with and get intuition of additive case

% vector of selection coefficients
%x_1 = [0.1; 0.3; 0.6; 0.1; 0.3; 0.6; 0.1; 0.3; 0.6; 0.1];
x_1 = [0.1; 0.3; 0.6];

m = size(x_1);
F = ones(m);

n = 10000;

%measure time for calculating it with correct math. notation
tic
for i = 1:n
    A_2 = x_1*F' + F*x_1';
end
toc
A_2

tic
for i = 1:n
    A_3 = x_1 + x_1';
end
toc
A_3


% matrix of fitness effects
A_1 = calc_A_additive(x_1);
A_1