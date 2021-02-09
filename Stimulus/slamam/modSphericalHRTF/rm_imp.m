function h = rm_imp(src,rec,walls,wtypes, ...
				sph_center,sph_rad,fs,c,taps,highpass);
%
% RM_IMP     calculate rectangular room source to receiver filter
%
%            RM_IMP(SRC,REC,WALLS,WTYPES, ...
%				SPH_CENTER,SPH_RAD,FS,C,TAPS,HIGHPASS) 
%	     returns a matrix of FIR taps which gives source to receiver
%            impulse response in a rectangular room. SRC is
%            SRC is a 1 by 3 matrix of the source location in
%            x-y-z space.  REC is a N by 3 matrix describing the
%            locations of N microphones.  WALLS is a 1 by 3 matrix of the
%	     [X Y Z] room dimensions.  The filter gain is normalized 
%	     for gain 1 distance unit from source.
%
%            WTYPES is a 1 by 6 matrix of surface material numbers.
%            See HELP ACOEFF for the material types (e.g. 26 = anec,
%	     27 = uniform (0.6) absorp coeff, 28 = uniform (0.2) absorp 
%	     coeff, and 29 = uniform (0.8) absorp coeff).  Specifically, 
%	     WTYPES(1) is the reflection coefficient for the wall
%            in the plane of X=0, WTYPE(2) for the wall at X=WALLS(1), 
%	     WTYPES(3) for the wall at Y=0, WTYPES(4) for the wall at
%	     Y=WALLS(4), etc.  
%	     NOTE: WTYPES can be a scalar, in which case all walls are
%	     assumed to be identical.
%
%	     SPH_CENTER and SPH_RAD are the center (in [x y z]) and radius of
%	     a rigid sphere placed in the room.  Note a scalar sphere center
%	     of 0 is automatically converted to [0 0 0] for simplicity.
%	     (default sph_center = [0 0 0] and sph_rad = 0 -- no sphere).
%
%            FS is the FIR filter sample rate (default = 10000 Hz).
%            C is the propagation velocity (default = 344.5).
%            TAPS is the number of taps to put in the impulse response.
%            (default = 2000).
%
%	     HIGHPASS is a flag that is 1 if the output impulse responses
%	     are to be highpass filtered to elimilate the DC component,
%	     otherwise no filtering is done.  Filter used is a two-zero,
%	     two-pole (z=1,1 and p=0.9889 +- j0.0110) butterworth.
%	     (default = 0, no highpass)
%
%	     RM_IMP uses the following MATLAB files:
%
%		IMP_GEN.M - Actually generates image source impulses
%				     and inserts into output impulse response.
%		SHAPEDFILTER.M - Used by IMP_GEN to create the
%				 non-integer impulse delay and to account
%				 for wall type and sphere diffraction.
%		SP.M	  \	
%		DELT.M	   \ - SP generates the sphere diffraction factor
%		HANK.M     /   using DELT, HANK and LEGPOLY
%		LEGPOLY.M /		

% Original by Mike O'Connell
% modified 9/96 by Jay Desloge


%
% Define variables as global to be used by IMP_GEN
%
 clear global
 global h src rec s_locations s_reflects 
 global sph_center sph_rad fs c taps ctap ctap2 fgains nfreq


%
% Check inputs and assign default values
%
if nargin < 3
	error('The source, receiver and wall locations must be specified');
end

if nargin < 10,
	highpass = 0,
end

if nargin < 9
	taps = 2000;	% Default output impulse response length
end

if nargin < 8
	c = 345;	% Default sound propagation velocity
end

if nargin < 7
	fs = 10000;	% Default sampling rate
	disp('room2src: assuming a 10000 Hz sampling rate.');
end

if nargin < 6		% Default setting = no sphere present
  	sph_center = [0 0 0];
	sph_rad = 0;
end
if sph_center == 0, sph_center = [0 0 0]; end

if nargin < 4
	wtypes = 26;    % Default wall type = all anechoic
	disp('room2src: assuming an anechoic enviornment');
end

if length(wtypes) == 1
	wtypes = wtypes * ones(1,6);
end

%
% find frequency dependent reflection coefficients for each wall
%
uniform_walls = 1;
for k = 1:6
	[alpha,freq] = acoeff(wtypes(k));  	% alpha = wall power absorption
						% freq = frequencies
	fgains(k,:) = sqrt(1-alpha);		% fgains = wall reflection 
        
	if abs(fgains(k,:)-fgains(k,1))>0, 	% If fgains depends on freq,
 		uniform_walls=0;		% set uniform_walls flag to 0.
	end
end

nfreq = freq/fs;				% freq as fraction of fs



%
% BEGIN CALCULATIONS
%

% 
% Part I: Initialization
%

% Set some useful values
 ctap=11;				% Center tap of lowpass to create
					% non-integer delay impulse
					% (as in Peterson)

 if (sph_rad==0)&(uniform_walls==1),	% Center tap of filter to account
  ctap2=1;				% for sphere presence and 
 else					% freq-dependent wall reflects.
  ctap2=17;				% if not present, omit.
 end
 num_rec=length(rec(:,1));		% number receivers


% Initialize output matrix (will later truncate to exactly taps length).
 h=zeros(taps+ctap+ctap2,num_rec);

%
% Part II: determine source image locations and corresponding impulse
% response contribution from each source.  To speed up process yet ease
% the computational burden, for every 10000 source images, break off and
% determine impulse response.
%

% The for determining source images is as follows.
%
%	1. Calculate maximum distance which provides relevant sources
%		(i.e., those that arrive within the imp_resp duration)
% 	2. By looping through the X dimension, generate images of
%		the (0,0,0) corner of the room, restricting the
%		distance below the presecribed level.
%	3. Use the coordinates of each (0,0,0) image to generate 8
%		source images
%	4. Generate corresponding number of reflections from each wall 
%	 	for each source image.

 dmax=ceil(taps*c/fs+max(walls));	% maximum source distance to be in
					% impulse response

 s_locations=ones(10000,3);		% Initialize locations and reflections
 s_reflects=ones(10000,6);		% matrices

 src_pts=[1 1 1;1 1 -1;1 -1 1;1 -1 -1;-1 1 1;-1 1 -1;-1 -1 1;-1 -1 -1].* ...
       (ones(8,1)*src);			% vector to get locations from
					% the (0,0,0) corner images

 Nx=ceil(dmax/(2*walls(1)));		% Appropriate number of (0,0,0)
					% images in either the +x of -x
					% directions to generate images
					% within dmax.

 loc_num=0;				% initialize location number index
 for nx=Nx:-1:0,			% loop through the images of (0,0,0)

  disp(sprintf('Stage %i',nx));

  if nx<Nx,			
   ny=ceil(sqrt(dmax*dmax - (nx*2*walls(1))^2)/(2*walls(2)));
   nz=ceil(sqrt(dmax*dmax - (nx*2*walls(1))^2)/(2*walls(3)));
  else					% Determine the number of y and z
   ny=0; nz=0;				% so that we contain dmax
  end

  X=nx*ones((2*ny+1)*(2*nz+1),1); 	% Form images of (0,0,0)
  Y=[-ny:ny]'*ones(1,2*nz+1); Y=Y(:); 
  Z=ones(2*ny+1,1)*[-nz:nz]; Z=Z(:); 
  if nx~=0,				% if nx~=0, do both +nx and -nx
   X=[-X;X]; Y=[Y;Y]; Z=[Z;Z];		% images of (0,0,0).
  end					
  Xw=2*walls(1)*X;
  Yw=2*walls(2)*Y;
  Zw=2*walls(3)*Z;

  for k=1:8,				% for each image of (0,0,0), get
   					% the eight source images and
					% number of relfects at each wall
   s_locs=zeros(length(X),3);
   s_refls=zeros(length(X),6);
   s_locs=[Xw,Yw,Zw]+ones(size(Xw))*src_pts(k,:);
   s_refls(:,1)=(src_pts(k,1)>0)*abs(X)+(src_pts(k,1)<0)*abs(X-1);
   s_refls(:,2)=abs(X);
   s_refls(:,3)=(src_pts(k,2)>0)*abs(Y)+(src_pts(k,2)<0)*abs(Y-1);
   s_refls(:,4)=abs(Y);
   s_refls(:,5)=(src_pts(k,3)>0)*abs(Z)+(src_pts(k,3)<0)*abs(Z-1);
   s_refls(:,6)=abs(Z);

   while (loc_num+length(s_locs(:,1)))>10000,	 % If the current source
    m=10000-loc_num;				 % image matrix has more than 
    s_locations((loc_num+[1:m]),:)=s_locs(1:m,:);% 10000 images, process to 
    s_reflects((loc_num+[1:m]),:)=s_refls(1:m,:);% get impulse response 
    done=imp_gen;			 	 % contributions.  

    loc_num=0;					 % Reset loc_num counter
    s_locs=s_locs(m+1:length(s_locs(:,1)),:);	 % and continue
    s_refls=s_refls(m+1:length(s_refls(:,1)),:);
   end

   s_locations((loc_num+[1:length(s_locs(:,1))]),:)=s_locs;
   s_reflects((loc_num+[1:length(s_locs(:,1))]),:)=s_refls;
   loc_num=loc_num+length(s_locs(:,1));		 % If current locations matrix
  end						 % has < 10000 images, keep
 end						 % building

 s_locations=s_locations(1:loc_num,:);		 % When all locations have 
 s_reflects=s_reflects(1:loc_num,:);		 % been generated, process
 done=imp_gen;				 % the final ones.


%
% Part III: Finalize output
%
 if highpass,
  [bhp,ahp] = butter(2,0.005,'high');	% Highpass filter, if desired
  for k=1:length(h(1,:)),
   h(:,k)=filter(bhp,ahp,h(:,k));
  end
 end
 h = h(1:taps,:);			% Restrict h to 'taps' length.
