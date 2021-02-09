function hout = shapedfilter(sdelay,freq,gain,th,rel,a,r,fs,ctap,ctap2);
%
% SHAPEDFILTER	shaped filter design for non-integer delay
%
%               SHAPEDFILTER(SDELAY,FREQ,GAIN,TH,REL,A,R,CTAP,CTAP2) is to
%		be used by IMPULSE_GENERATE to generate the
%		(2*CTAP-1)+(2*CTAP2-1)-1 long impulse responses (in each 
%		row) that produce the non-integer delay (SDELAY) and 
%		freq-dependent wall-reflection/sphere diffraction effects 
%		for each source image location.
%
%		The inputs are:
%		
%		SDELAY a vector of non-integer delays for each source image.
%
%		FREQ, GAIN are vectors specifying the frequency dependent 
%		wall reflectin effects for each source image.
%
%		TH is a vector of source-to-receiver locations for each image.
%
%		REL is a vector of the ratio of image-source magnitude to
%		direct path magnitude.
%
%		A, R are the sphere radius and the distance from the sphere 
%		center to the receiver, respectively.
%
%		FS is the sampling frequency.
%		CTAP is the center tap of the filter for non-integer delay
%		CTAP2 is the center tap of filter to account for the
%		      frequency-dependent wall reflections and 
%		      sphere diffraction. 

% Originally by Mike O'Connell
% Revised by Jay Desloge 9/19/96


%
% Parse inputs
%
 if nargin < 6,
        ctap2 = 17;
 end

 if nargin < 5
	ctap = 11;
 end

 if nargin < 4
	error('Frequencies, gains and sdelay must be specified');
 end

 if sdelay < 0
	error('The sample delay must be positive');
 end

%
% design the non-integer delay filter
%
 ntaps = 2*ctap-1;
 N=ctap-1;
 fc=0.9; 
 h=0.5*fc*(1+ ...
    cos(pi*(ones(size(sdelay))*[-N:N]-sdelay*[0,ones(1,ntaps-1)])/N)).* ...
    sinc(fc*(ones(size(sdelay))*[-N:N]-sdelay*[0,ones(1,ntaps-1)]));

%
% make sure that freq are is vectors
%
 freq=freq(:).';


%
% design and incorporate the wall filter/sphere filter if it is needed 
% (i.e., is ctap2>1).  Otherwise, just scale impulse resp. appropriately.
%

 if ctap2>1,
   df = [0:ctap2-1]*(pi/(ctap2-1)); 		% Determine FFT points

   freq = [ -eps 2*pi*freq pi];	% frequencies at which gain is defined
   gain = [gain(:,1) gain gain(:,length(gain(1,:)))];
   G = interp1(freq.',gain.',df).';
				% Interpolate reflection frequency-dependence
				% to get gains at FFT points.


  % If the sphere is present, incorporate it
  % (done for all location simultaneously).
  %
   if a > 0,
    k=df*fs/345;			% wave number
    direct=zeros(length(th),length(k));% direct-wave component
    yc=zeros(size(direct));		% direct+scattered

    T=sp(k(2:length(k)),[0:45]*pi/45,a,r);
					% Look-up table for weaker reflections
     
    if sum(find(rel<0.01))>0,    %if there are weaker reflections             
    direct(rel<0.01,:) = [exp(-j*r*(-cos(th(rel<0.01)))*k)];
    yc(rel<0.01,2:length(k))=T(1+round(45*th(rel<0.01)/pi),:);
					% For weaker reflections, find nearest
					% direct and direct+scattered
    end
    
    if sum(rel>=0.01)>0,
     direct(rel>=0.01,:) = [exp(-j*r*(-cos(th(rel>=0.01)))*k)];
     yc(rel>=0.01,2:length(k)) = sp(k(2:length(k)),th(rel>=0.01),a,r);
					% For stronger reflections, determine
					% direct and direct+scattered from sum.
    end

    yc(:,1)=ones(length(th),1);		% assign 1 to DC scattered component 

    if(size(G,1) ~= size(yc,1))
        G = G';
    end
    G= yc.*G./direct;			% modify the gf above (which represents
					% the direct-wave component) by 
					% (direct+scattered./direct) to 
					% convert into sphere-scattered filter
   end

  %
  % Combine the non-integer delay filter and the wall/sphere filter
  %
   G(:,ctap2)=real(G(:,ctap2));
   G=[G,fliplr(conj(G(:,2:ctap2-1)))];    	% Transform into appropriate
					% wall into transfer function.
   gt=real(ifft(G.'));			% IFFT to get imp-resp
   
   g=[0.5*gt(ctap2,:);gt(ctap2+1:2*ctap2-2,:); ...
      gt(1:ctap2-1,:);0.5*gt(ctap2,:); ...
      zeros(2*ctap-2,length(gt(1,:)))];
   G=fft(g);				% Zero-pad and FFT

   H=fft([h.';zeros(2*ctap2-2,length(h(:,1)))]);
					% Zero-pad and FFT delay filter

   HOUT=H.*G;				% Convolve filters
   hout=real(ifft(HOUT)).';		% Obtain total impulse response
 
 else
 
   hout=h.*(gain(:,1)*ones(size(h(1,:))));
					% Scale impulse response only if
					% Wall reflects are freq-indep and
					% if sphere not present.
 end