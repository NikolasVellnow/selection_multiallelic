% Script to play around with and get intuition of nfds

% vector of selection coefficients
%x_1 = [0.1; 0.3; 0.6; 0.1; 0.3; 0.6; 0.1; 0.3; 0.6; 0.1];
x_1 = [0.01; 0.03; 0.06];

m = size(x_1);
F = ones(m);

n = 10000;

F*F'

%measure time for calculating it with correct math. notation
tic
for i = 1:n
    A_2 = (F + x_1) * (F + x_1)' - F*F';
end
toc
A_2

tic
for i = 1:n
    A_3 = x_1 .* x_1';
end
toc
A_3


% matrix of fitness effects
A_1 = calc_A_multiplicative(x_1);
A_1