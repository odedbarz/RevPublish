function [C, P] = corrcoef_array(X, Y)
%
%   function [C, P] = corrcoef_array(X, Y)
%
% Description:
% Pearson correlation coefficient over each column in the matrices X & Y.

nt = size(X,2);
C = zeros(nt, 1);
P = zeros(nt, 1);

% Avoid NaN correlations when correlating between a vector(s) of zeros
eps = 1e-12;     

for k = 1:nt
    [ck, pk] = corrcoef(X(:,k)+eps, Y(:,k)+eps);
    C(k) = ck(1,2);
    P(k) = pk(1,2);    
end