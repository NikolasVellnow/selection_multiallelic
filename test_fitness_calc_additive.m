% Script to play around with and get intuition of additive case

% vector of selection coefficients
%s_1 = [0.1; 0.3; 0.6; 0.1; 0.3; 0.6; 0.1; 0.3; 0.6; 0.1];
s_1 = [0.1; 0.3; 0.6];

m = size(s_1);
F = ones(m);

n = 10000;

%measure time for calculating it with correct math. notation
tic
for i = 1:n
    A_2 = s_1*F' + F*s_1';
end
toc
A_2

tic
for i = 1:n
    A_3 = s_1 + s_1';
end
toc
A_3


% matrix of fitness effects
A_1 = calc_A_additive(s_1);
A_1