function [A] = calc_A_heterozygote(sigma, n)
% Calculates entries of matrix of fitness effects A_ij(x) with heterozygote
% advantage sigma common to all heterozygotes and number of n alleles

F = ones([n,1]);
I = eye(n);


% matrix of fitness effects
A = sigma * (F*F' - I);

end