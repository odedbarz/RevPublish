function leg = legpoly(n,x)
%
% LEGPOLY	Legrendre polynomial
%
%       LEGPOLY(n,X) calculates the first n terms (0:n-1) of the
%       Legrendre polynomial of the values of vector X.
%       The result is a matrix with n
%       rows and length(X) columns.
%
% modified from original function leg from pmz
% mpo 6/7/95

x=x(:).';

leg(1,:) = ones(1,length(x));
leg(2,:) = x;

for k = 3:n
         leg(k,:) = ( x.*(2*k-3).*leg(k-1,:) - (k-2).*leg(k-2,:) ) / (k-1);
end

