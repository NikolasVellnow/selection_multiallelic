% Script to play around with and get intuition of heterozygote advantage

% vector of selection coefficients
%x_1 = [0.1; 0.3; 0.6; 0.1; 0.3; 0.6; 0.1; 0.3; 0.6; 0.1];
x_1 = [0.01; 0.03; 0.06];

sigma = 0.03;

n = length(x_1);
F = ones([n,1]);
I = eye(n);

m = 10000;

(F*F' - I)



%measure time for calculating it with correct math. notation
tic
for i = 1:m
    A_2 = sigma * (F*F' - I);
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
A_1 = calc_A_heterozygote(sigma, n);
A_1