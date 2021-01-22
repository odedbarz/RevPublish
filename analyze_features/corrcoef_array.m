function [C, P] = corrcoef_array(X, Y)
%
%   function [C, P] = corrcoef_array(X, Y)
%
% Description:
% Pearson correlation coefficient over each column in the matrices X & Y.

nt = size(X,2);
C = zeros(nt, 1);
P = zeros(nt, 1);

for k = 1:nt
    [ck, pk] = corrcoef(X(:,k), Y(:,k));
    C(k) = ck(1,2);
    P(k) = pk(1,2);
    
end