function Xw = zca(X, d)
%
% function Xw = zca(X)
%
% Description:
% Whitening (sphering) of each column in X. That means that each column in
% X has zero mean and a variance of 1.0.

if 2 > nargin 
    d = 1;
end

Xw = (X-mean(X, d))./(eps + std(X, [], d));