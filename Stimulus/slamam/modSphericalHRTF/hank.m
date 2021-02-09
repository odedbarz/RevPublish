function hank = hank(n,kr)
%
% HANK
%
%          HANK(N,KR) calculates the first N
%          (0:n-1) Hankel functions with argument KR.
%

% modified from the function hd.m from pmz
% mpo 6/6/95
% jgd 9/18/96 

 % 
 % Initialize the series
 %
  kr=kr(:)';
  h0 = cos(kr)./kr + (sin(kr)./kr)*j;		% Hankel Function
  hank(1,:) = sin(kr)./kr - (cos(kr)./kr)*j;
  hank(2,:) = hank(1,:)./kr - h0;

 %
 % Iterate to generate hankel function.  As long as index, k, is less than
 % or equal to the maximum summation index for input KR(i), keep generating
 % hankel functions.  Otherwise, truncate hankel functions to 0 (as then
 % are not used in the sum.
 %
 for k=3:max(n);
  % Handle single-frequency case (KR has one element) differently to
  % allow for vectorizing of operations.
   if length(n)>1,
    hank(k,n>=k) = (2*k-3)*hank(k-1,n>=k)./kr(n>=k) - hank(k-2,n>=k);
    hank(k,n<k) = zeros(1,sum(n<k));
   else
    hank(k) = (2*k-3)*hank(k-1)./kr - hank(k-2);
   end   
 end

