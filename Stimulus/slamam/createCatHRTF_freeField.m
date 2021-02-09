%
% createCatHRTF_test.m
%
% computes one HRTF using spherical/image model
%
% Goal: adjust cat's head radius, ear position etc... 
%  to acheieve ideal ITD-azimuth function for cat
%
%

addpath c:\users\devores\catHRTF\modSphericalHRTF

walls = [11 13 3]; %(in meters)

%wtypes = [20 20 20 20 20 20];
wtypes =  [7 7 7 7 14 16]; %reverberant wall types 

%create a room with no sphere
sph_cent = [4.7 5.4 1.4]; %sphere is in center of room
sph_r = 0;

%one receiver at the center of the room
recs(1,:) = sph_cent;

f_samp = 50000;
c_snd = 344.5    

num_taps = 0 %150000; %let model figure out how many taps
highpass = 0; %don't filter out DC component
dsply = 1; %show progress
err = 0; %no error introduced

%matrix of sources 
angle = [0];
distance = [ 0.5 1 1.5 3];

for ds=1:3
    [xt,yt,zt] = polar2rect(distance(ds),angle(1));
    sources(ds,:) = sph_cent + [xt yt zt];
end

src_list = {'0d_50cm','0d_100cm', '0d_150cm', '0d_300cm'};


for i= 1:1 %length(sources)
    disp(sprintf('Calculating impulse response for source %s',char(src_list(i))));
    [H_OUT LEAD_ZEROS] = room_impulse(sources(i,:),recs,walls,wtypes,sph_cent,sph_r,f_samp,c_snd,num_taps,highpass,dsply,err);
    eval(sprintf('save freeField_%s *',char(src_list(i))));
    clear H_OUT LEAD_ZEROS
end
