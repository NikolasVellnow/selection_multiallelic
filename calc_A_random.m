function [A] = calc_A_random(n, s0)
%Calculates entries of matrix of symmetric random fintess effects
%based on the number of alleles and the strength of selection

% matrix of random fitness effects
B = randn(n);
A = s0 * (B + B');

end