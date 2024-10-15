% Script to play around with and get intuition of heterozygote advantage

% vector of selection coefficients
%x_1 = [0.1; 0.3; 0.6; 0.1; 0.3; 0.6; 0.1; 0.3; 0.6; 0.1];
x_1 = [0.01; 0.03; 0.06];

lambda = 0.5;

n = length(x_1);


m = 10000;

Z = randn(n);

Z



%measure time for calculating it with correct math. notation
tic
for i = 1:m
    A_2 = sqrt(lambda)*(Z + Z')/2;
end
toc
A_2

% is there a faster way?
%tic
%for i = 1:m
    %A_3 = sigma * (F*F' - I);
%end
%toc
%A_3


% matrix of fitness effects
A_1 = calc_A_fluctuating(lambda, n);
A_1