function [A] = calc_A_fluctuating(lambda, n)
% Calculates entries of matrix of fitness effects A_ij(x) with heterozygote
% advantage sigma common to all heterozygotes and number of n alleles

Z = randn(n);

% matrix of fitness effects
A = sqrt(lambda)*(Z + Z')/2;
end