function yc = sp(k,theta,a,r)
%
% SP	complex pressure just off the surface of a sphere
%
%    SP(k,theta,a,r) calculates the complex pressure factor
%    on the surface sphere.  In this program, 
%          
%           k = N by 1 vector of wave_number = omega/c (c=speed of sound, m/s) 
%           theta = M by 1 vector angles between the wave's direction 
%                   of arrival and point of measurement 
%           a = sphere radius (m)
%           r = distance from center of sphere to point of measurement (m)
%
%    Output is a M by N matrix where each row is the complex pressure of
%    a given angle as a function of frequency.
%

%  modified function cp.m from pmz; spherepressure.m from mpo
%  jgd 9/18/96

%
% Initialize output signal.  Assign output to 1 for all k=0 (DC) values,
% and perform iteration on non-zero k values.
%
 yc = ones(length(theta),length(k));
 k = k(:);
 i_non_zero = find(k~=0);

%
% Form vectors of k*a and k*r
%
 ka=k(i_non_zero)*a;
 kr=k(i_non_zero)*r;

%
% Evaluate Delta Parameter, Hankel function and the LeGendre Polynomials
% at the desired frequencies, angles, and summation indices used in
% Morse Equation 29.9
%
 [dm,n]=delt(ka);		% Determine delta parameter for Morse 29.9
				% and summation index for each frequency.
				% dm is a matrix, with row = summation index
				% and col = freq.

 hnk=hank(n,kr);		% Determine hankel function at desired freqs
				% for summation indices determined above.
				% hnk is a matrix with row = summation index
				% and col = freq.

 h_terms=hnk.*exp(-j*dm).*sin(dm);	% Combine delta terms and hankel 
					% functions into the terms used
					% by Morse 29.9. 

 leg = legpoly(max(n),-cos(theta)); 	% Generate Legendre polynomials, 
					% Lm(cos(theta)), at the max summation
					% indices and for the desired angles.
					% leg is a matrix with row = summation
					% index and col = theta.
 %
 % Generate Direct wave
 %
  s = exp(j*(-cos(theta(:)))*kr(:).');

 %
 % Loop through the frequencies.  At each frequency, generate the scattered
 % wave using Equation 29.9 and add to the direct wave.
 %
  for freq=1:length(k(i_non_zero)),
   m = [0:n(freq)-1]';
   x1 = ((-(2*m+1).*(i.^(m+1)).*h_terms(1:n(freq),freq))* ...
		ones(1,length(theta))).*leg(1:n(freq),:);
   x2 = ones(1,n(freq))*x1;

   % Add scattered term to direct term to get output.
   s(:,freq) = s(:,freq) + x2.';
  end

 %
 % Assign output pressure to non-zero components of yc.
 %
  yc(:,i_non_zero) = conj(s);







