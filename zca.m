function Xw = zca(X)
%
% function Xw = zca(X)
%
% Description:
% Whitening (sphering) of each column in X. That means that each column in
% X has zero mean and a variance of 1.0.

Xw = (X-mean(X))./(eps + std(X));