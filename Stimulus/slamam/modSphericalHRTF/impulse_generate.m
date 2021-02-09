function done=impulse_generate
%
% IMPULSE_GENERATE
%
%	IMPULSE_GENERATE is called by the function ROOM_IMPULSE.
%	It take a set of source images (S_LOCATIONS) and generates
%	the corresponding impulse response component between each image
%	and a set of receivers (REC) and inserts this into the overall
%	impulse response for each receiver (H).
%
%	All variables are passed as global variables from ROOM_IMPULSE.
%	Other variables passed are:
%	
%	S_REFLECT = the number of reflections at each wall ecperienced by
%	each source image.
%
%	FGAINS, NFREQ = the wall reflection freq-dependent gains and
%	corresponding frequencies.
%
%	SPH_CENTER, SPH_RAD = center and radius of a rigid sphere,
%	diffraction about which affects the source-image to microphone
%	impulse response.
%	
%	FS = sampling frequency.
%	C = sound propoagation velocity.
%	TAPS = desired length of output impulse response.
%	CTAP = center tap of filter to approximate the non-integer sample
%	       delay when travelling from source to receiver.
%	CTAP2 = center tap of filter to approximate the combined effects of
%		freq-dependent wall reflections and diffraction about sphere.
%	LEAD_ZEROS = the number of leading zeros to strip from the impulse 
%		response.
%

% J. Desloge 9/19/96

% 
% Global variables from ROOM_IMPULSE
%
 global h src rec s_locations s_reflects fgains nfreq
 global sph_center sph_rad fs c taps ctap ctap2 lead_z dspy

%
% Part I: Form variables to be used in impulse response generation
%

% Determine the overall source gains (based on number of reflections
% through each wall) for each source location.
 gains = ones(length(s_locations(:,1)),length(nfreq));		
 for wall = 1:6
  gains = gains .* ...
          ((ones(size(s_locations(:,1)))*fgains(wall,:)).^ ...
             (s_reflects(:,wall)*ones(size(nfreq)))); 
 end	
 ind = find(sum(gains.').'>0);
 if isempty(ind), ind = 1; end
 s_locations = s_locations(ind,:);
 s_reflects = s_reflects(ind,:);
 gains = gains(ind,:);

% If the sphere is actually present, calculate the distances between the 
% sources and the sphere center and between the receivers and the sphere 
% center -- to be used to determine the angle between source and receiver
 num_rec=length(rec(:,1));
 s_dist_sph = zeros(length(s_locations),1);
 r_dist_sph = zeros(num_rec,1);
 if sph_rad>0,
  dist = s_locations-(ones(size(s_locations(:,1)))*sph_center);
  s_dist_sph=((dist.*dist)*[1;1;1]).^(0.5);	% Source distance

  dist = rec-(ones(num_rec,1)*sph_center);
  r_dist_sph=((dist.*dist)*[1;1;1]).^(0.5); 	% Receiver distance
  if (min(r_dist_sph)<(sph_rad-eps)), 
   error('Receiver located inside sphere!!!!'); 
  end;
 end

%
% Part II: Determine impulse response components for each receiver 
% and insert into overall impulse response

 for kk = 1:num_rec, 
   if dspy,
		disp(sprintf('Receiver %i out of %i',kk,num_rec));
	end

  % If sphere present, calculate angle of arrival between the 
  % source-sphere_center ray and the receiver-sphere_center array 
  % (for use in figuring effect of sphere).
   th=zeros(length(s_locations(:,1)),1);
   if sph_rad>0,
    u_source=(s_locations - ones(size(s_locations(:,1)))*sph_center)./ ...
                (s_dist_sph*[1 1 1]); 	% unit vector sph_center -> source
    u_rec=(rec(kk,:)-sph_center)./(r_dist_sph(kk)*[1 1 1]);
  					% unit sph_center -> rec(kk)
    th=acos(u_source*u_rec.');		% Angle between source/receiver
   end

  % Calculate the distance between the source locations and the
  % current receiver location.
   sr_dist = zeros(length(s_locations),1);
   dist = s_locations-(ones(size(s_locations(:,1)))*rec(kk,:));
   sr_dist=((dist.*dist)*[1;1;1]).^(0.5);

  % Calculate samples delay for source to current receiver
   thit = ctap + ctap2 - lead_z + ...
				sr_dist/c*fs; 			% number of samples delay to travel to rec
   ihit = fix(thit);					% round down	
   
   
   %!!!!!!modified by sasha 2-3-05!!!!!
   %forces all non-integer delays to be zero so that we don't use the
   %shaped filter
   %fhit = thit - ihit;				% non-integer extra samples
   fhit = zeros(size(ihit));
    
    
  % eliminate locations that are too far away to enter into impulse response
   v=(ihit<=taps+ctap+ctap2);
   

  % if there are actual impulses to calculate, do it.  Otherwise, skip.
  
   if sum(v)>0,
  
  % Calculate impulse strength re. direct path (for purposes of simpifying
  % sphere calculations).
    direct_dist=norm(src-rec(kk,:),2);
    rel_direct=direct_dist*(max(gains.').')./sr_dist;
  

  % Initialize temporary impulse response vector
    ht = zeros(taps+2*ctap+1+2*ctap2+1,1);
   
    ht_ind = (ihit(v)*ones(1,2*ctap-1+2*ctap2-1-1))+ ...
		(ones(size(ihit(v)))*[-ctap-ctap2+1+1:ctap+ctap2-1-1]);
				% Indices into ht.  Each row corresonds to 
				% one source image location, with the center
				% determined by ihit. Within a row, there are
				% (2*ctap-1)+(2*ctap2-1)-1 values
				% that account for non-integer delay, fhit and
				% for the freq-dep wall reflections/sphere
				% diffraction.

  % For each source location, determine the impulse response.
    h_temp = 1./(sr_dist(v)*ones(1,2*ctap-1+2*ctap2-1-1)) ...
        .*(shapedfilter(fhit(v),nfreq,gains(v,:),th(v),rel_direct(v), ...
			   sph_rad,r_dist_sph(kk),fs,ctap,ctap2));
				% form filter to incorporate frequency gains,
				% non-integer delay and scattering off of 
				% rigid sphere.

  % Add the impules response segments into the overall impulse response.
  %keyboard
    for k=1:sum(v),
        if any(ht_ind(k,:)<1), keyboard, end
        ht(ht_ind(k,:),1)=ht(ht_ind(k,:),1)+h_temp(k,:)';
    end

   % Add into overall impulse response matrix
    h(:,kk)=h(:,kk)+ht(1:taps+ctap+ctap2); 
   end
 end
done = 1;






