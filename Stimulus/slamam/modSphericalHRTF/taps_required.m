function [num_taps] = taps_required(srcs,mics,walls,wtypes,f_samp,c_snd);
%
% [num_taps] = taps_required(srcs,mics,walls,wtypes,f_samp,c_snd);
%
% This function determines the 

%
% find frequency dependent reqlection coefficients for each wall and
% choose the maximum one
%
if length(wtypes) == 1
	wtypes = wtypes * ones(1,6);
end
wgains = zeros(6,1);
for k = 1:6,
	[alpha,freq] = acoeff(wtypes(k));  	% alpha = wall power absorption
						                % freq = frequencies
	wrefl(k) = max((1-alpha));	        % wgains = max wall reflection 
end

% max distance mic is reference.
dists = sqrt(((mics-ones(size(mics(:,1)))*srcs).^2)*ones(3,1));
[max_dists,i_dist]=max(dists);
m_ref = mics(i_dist,:);

%
% Generate scale factor direct wave.
% For each dimension (X,Y,Z), step out until scale factor (based
% on distance and wall reflection coeffs) is 1/1000000 of ref.
%
dir = srcs-m_ref;
ref_scale = 1./sum(dir.^2);
max_dist = zeros(3,1);
for dim = 1:3,
    inds = (dim-1)*2 + [1:2];
    step = zeros(1,3);
    step(dim) = 2*walls(dim);
    refl = prod(wrefl(inds));
    
    k = 1;
    scale = (refl^k)./sum((dir+k*step).^2);
    while scale>0.000001*ref_scale,
        k = k+1;
        scale = (refl^k)./sum((dir+k*step).^2);
    end        
    k = k+1;    % give an extra step just for good measure.
    max_dist(dim) = norm(dir+k*step);
end

num_taps = ceil(f_samp*(max(max_dist)/c_snd));
num_taps = ceil(1.5*num_taps);
