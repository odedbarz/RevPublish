%
% createCatHRTF_nosphere_2recs.m
%
% computes one HRTF using spherical/image model
%
% Goal: creates 2-receiver bilateral room impulse responses using a
% modified version of the room-image method that does not allow for
% non-integer delays
%

%addpath u:\delgutte\devores\catHRTF\modSphericalHRTF
addpath .\modSphericalHRTF     % Oded, 06/12/2019

walls = [11 13 3]; %(in meters)

%wtypes = [20 20 20 20 20 20];
wtypes =  [7 7 7 7 14 16]; %reverberant wall types 

%create a room with no sphere
sph_cent = [4.7 5.4 1.4]; %sphere is in center of room
sph_r = 0;

%two receivers 10 centimeters apart
recs(1,:) = sph_cent + [-0.06 0 0];
recs(2,:) = sph_cent + [0.06 0 0];


f_samp = 50000;
c_snd = 344.5    

num_taps = 0 %150000; %let model figure out how many taps
highpass = 0; %don't filter out DC component
dsply = 1; %show progress
err = 0; %no error introduced

%matrix of sources 
angle = [-3.8] % 15 30 45 60 75 90 -15 -30 -45 -60 -75 -90];
%angle = [];
%distance = [1 1.5 3 6];
distance = [0.5 1 1.5 3];

for ds=1:4 %4
    for ag=1:1 %13
        [xt,yt,zt] = polar2rect(distance(ds),angle(ag));
        (ds-1)*1 + ag
        sources((ds-1)*1 + ag,:) = sph_cent + [xt yt zt];
    end
end

src_list = {'-4d_50cm','-4d_100cm','-4d_150cm','-4d_300cm'};
        
%src_list = {'-4d_50cm' '15d_50cm' '30d_50cm' '45d_50cm' '60d_50cm'  '75d_50cm' '90d_50cm' '-15d_50cm' '-30d_50cm' '-45d_50cm' '-60d_50cm' '-75d_50cm' '-90d_50cm','-4d_100cm'}; %...
%        '0d_150cm' '15d_150cm' '30d_150cm' '45d_150cm' '60d_150cm'  '75d_150cm' '90d_150cm' '-15d_150cm' '-30d_150cm' '-45d_150cm' '-60d_150cm' '-75d_150cm' '-90d_150cm' ...
%        '0d_300cm'  '15d_300cm' '30d_300cm' '45d_300cm' '60d_300cm' '75d_300cm' '90d_300cm'  '-15d_300cm' '-30d_300cm' '-45d_300cm' '-60d_300cm' '-75d_300cm' '-90d_300cm' ...
%    '0d_600cm' '15d_600cm' '30d_600cm' '45d_600cm' '60d_600cm'  '75d_600cm' '90d_600cm' '-15d_600cm' '-30d_600cm' '-45d_600cm' '-60d_600cm' '-75d_600cm' '-90d_600cm'};


for i= 1:4 %length(sources)%1:length(sources)
    disp(sprintf('Calculating impulse response for source %s',char(src_list(i))));
    [H_OUT LEAD_ZEROS] = room_impulse(sources(i,:),recs,walls,wtypes,sph_cent,sph_r,f_samp,c_snd,num_taps,highpass,dsply,err);
    
    % Oded, 06/12/2019: skip the saving 
    %eval(sprintf('save hrtf_itd_%s *',char(src_list(i))));
    
    clear H_OUT LEAD_ZEROS
end
