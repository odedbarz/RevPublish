function dcor = distcorr(x,y)
%
% This function calculates the distance correlation between x and y.
% Reference: http://en.wikipedia.org/wiki/Distance_correlation 
% Date: 18 Jan, 2013
% Author: Shen Liu (shen.liu@hotmail.com.au)
%
% Changes:
% * I have changed the doubly-centered-distance matrix procedure into one inner
%   function;
%
% * I have changed the doubly-centered-distance matrix procedure to use MATLAB's 
%   vectorization and so to avoid creating huge matrices (OBZ).
%
% * Each COLUMN is one sample 

import distcorr.*

% The two random variables x & y must have the same number of samples (columns)
assert( size(x,2) == size(y,2),...
    'Error @ distcorr.distcorr(x,y): Inputs must have the same number of samples (columns)!');

% Delete rows containing unobserved values
N = any([isnan(x) isnan(y)],2);
x(N,:) = [];
y(N,:) = [];

% Calculate doubly centered distance matrices for x and y
A = doubly_centered_dist(x);
B = doubly_centered_dist(y);

% Calculate empirical distance *covariance*   
%dcov = sum(sum(A.*B))/(size(mrow,1)^2);
dcov = distcov(A, B);

% Calculate empirical distance *variance* for X & Y    
% dvarx = sum(sum(A.*A))/(size(mrow,1)^2);
% dvarx = max(0, dvarx);  % the variance is always semi-positive
dvarx = distcov(A, A);

% dvary = sum(sum(B.*B))/(size(mrow,1)^2);
% dvary = max(0, dvary);  % the variance is always semi-positive
dvary = distcov(B, B);

% Calculate the empirical distance *correlation*    
if 0 == dvarx * dvary
    dcor = 0;   % by definition of dCor
else
    dcor = dcov/sqrt(dvarx*dvary);
end





