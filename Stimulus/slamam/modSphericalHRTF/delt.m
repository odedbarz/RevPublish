function [delta,n] = delt(ka)
%
% DELT
%
%          [DELTA,N]=DELT(KA).  For each input value, KA(i), DELT 
%	   calculated the 'delta' parameter values at N(i) summation
%	   indices.  N(i) is the maximum summation index for the
%	   ith term and is determined as the index such that the delta
%	   value corresponding to N(i)+1 is below 1e-14.
%

% modified from the function hd.m from pmz
% mpo 6/6/95
% jgd 9/18/96 

% 
% Initialize the series
%
 ka=ka(:)'; n=zeros(size(ka));

 h0 = cos(ka)./ka + (sin(ka)./ka)*j;		% Hankel Function
 hank(1,:) = sin(ka)./ka - (cos(ka)./ka)*j;
 hank(2,:) = hank(1,:)./ka - h0;

 dhank(1,:) = h0 - hank(1,:)./ka;		% Hankel Function derivative
 dhank(2,:) = hank(1,:) - hank(2,:)*2./ka;	

 delta(1:2,:)=angle(dhank/j);			% Delta term

%
% Iterate through summation index, indx, to generate hankel functions, 
% derivatives and delta terms for all KA(i) at the current index.  
% When the series of delta's corresponding to KA(i) is 
% small enough (less than 1e-14, in the noise floor), truncate to zero and 
% assign the corresponding maximum index N(i).  When all delta series for 
% all KA(i) are below the threshold, break out of the loop. 
%
 done=0; indx=3;
 while ~done
  %
  % Hankel functions and derivatives.
  %
   hank(indx,:) = (2*indx-3)*hank(indx-1,:)./ka - hank(indx-2,:);
   dhank(indx,:) = hank(indx-1,:) - hank(indx,:)*indx./ka;

  %
  % Delta terms.
  %
   delta(indx,:)=angle(dhank(indx,:)/j);

  %
  % Assign maximum summation index to the delta series that drop
  % below the threshold (1e-14) for the first time.
  %
   n(find((n==0)&(abs(delta(indx,:))<1e-14)))=indx* ...
			ones(1,sum((n==0).*(abs(delta(indx,:))<1e-14)));

  %
  % Truncate any delta's below threshold to zero.
  %
   delta(indx,:)=(abs(delta(indx,:))>=1e-14).*delta(indx,:);

  %
  % If all delta series [corresponding to all inputs KA(i)] have fallen
  % below threshold,  exit while loop.
  %
   if max(abs(delta(indx,:)))==0, done=1; end
   indx=indx+1;
 end

