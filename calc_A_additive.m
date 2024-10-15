function [A] = calc_A_additive(s)
%Calculates entries of matrix of additive fitness effects A_ij(x)
%based on selection coefficients in column vector s

% matrix of additive fitness effects
A = s + s';     % corresponds to s*F' * F*s' but is 3x faster

end