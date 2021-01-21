function [C, P] = corrcoef_array(X, Y)

xdim = size(X,2);
C = zeros(xdim, 1);
P = zeros(xdim, 1);

for k = 1:xdim
    [ck, pk] = corrcoef(X(:,k), Y(:,k));
    C(k) = ck(1,2);
    P(k) = pk(1,2);
    
end