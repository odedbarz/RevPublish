function dC = distcov(A, B)
%
% Calculates the empirical distance covariance
% 

dC = 1/size(A,1) * sqrt(sum( A(:) .* B(:) ));

% The variance is always semi-positive
dC = max(0, dC);

