function [A] = calc_A_nfds(x,c)
%Calculates entries of matrix of fitness effects A_ij(x)
%based on allele frequencies in column vector x and proportionality constant c

% matrix of genotype frequencies
G = x * x';

% matrix of negative frequency-dependent fitness effects
A = (-c) * G;
end