%
% createCatHRTF_nosphere_2recs.m
%
% computes one HRTF using spherical/image model
%
% Goal: creates 2-receiver bilateral room impulse responses using a
% modified version of the room-image method that does not allow for
% non-integer delays
%

addpath 'C:\Documents and Settings\slamam\My Documents\My Research\IC Research\HRTFs\catHRTF\modSphericalHRTF'
 
walls = [11 13 3]; %(in meters)

%wtypes = [26 26 26 26 26 26];
wtypes =  [7 7 7 7 14 16]; %reverberant wall types 

%create a room with no sphere
sph_cent = [4.7 5.4 1.4]; %sphere is in center of room
%sph_r = 0.05999; that's for rabbits and cats
sph_r = 0.10999;    % for humans 22 cm interaural distance

%two receivers 22 centimeters apart
% recs(1,:) = sph_cent + [-0.06 0 0]; for bunnies and cats
% recs(2,:) = sph_cent + [0.06 0 0];
recs(1,:) = sph_cent + [-0.11 0 0];
recs(2,:) = sph_cent + [0.11 0 0];

f_samp = 50000;
c_snd = 344.5;    

num_taps = 0 %150000; %let model figure out how many taps
highpass = 0; %don't filter out DC component
dsply = 1; %show progress
err = 0; %no error introduced

%matrix of sources 
angle = 0;%-90:5:90;
distance = [1.5 3];%[0.5 1 1.5 3];
    
for ds=1:length(distance) %4
    for ag=1:length(angle)
        [xt,yt,zt] = polar2rect(distance(ds),angle(ag));
        (ds-1)*length(angle) + ag
        sources((ds-1)*length(angle) + ag,:) = sph_cent + [xt yt zt];
    end
end

%src_list = {'-0d_50cm' '15d_50cm' '30d_50cm' '45d_50cm' '60d_50cm'  '75d_50cm' '90d_50cm' '-15d_50cm' '-30d_50cm' '-45d_50cm' '-60d_50cm' '-75d_50cm' '-90d_50cm', ...
%        '0d_100cm' '15d_100cm' '30d_100cm' '45d_100cm' '60d_100cm'  '75d_100cm' '90d_100cm' '-15d_100cm' '-30d_100cm' '-45d_100cm' '-60d_100cm' '-75d_100cm' '-90d_100cm', ...
%        '0d_150cm' '15d_150cm' '30d_150cm' '45d_150cm' '60d_150cm'  '75d_150cm' '90d_150cm' '-15d_150cm' '-30d_150cm' '-45d_150cm' '-60d_150cm' '-75d_150cm' '-90d_150cm', ...
%        '0d_300cm'  '15d_300cm' '30d_300cm' '45d_300cm' '60d_300cm' '75d_300cm' '90d_300cm'  '-15d_300cm' '-30d_300cm' '-45d_300cm' '-60d_300cm' '-75d_300cm' '-90d_300cm'};

 
figure(100);
for ds=1:length(distance) %4
for ag=1:length(angle)
[xt,yt,zt] = polar2rect(distance(ds),angle(ag));
scatter3(xt+sph_cent(1),yt+sph_cent(2),zt+sph_cent(3),'bo');
hold on
end
end
scatter3(recs(1,1), recs(1,2), recs(1,3),'r.','LineWidth',2)
scatter3(recs(2,1), recs(2,2), recs(2,3),'r.','LineWidth',2)
axis([0 11 0 13 0 3]);
disp('Please verify source arrngements in plot and then hit ENTER');
pause;

i=0;
for ds=1:length(distance) 
    for ag=1:length(angle)  
    i=i+1;    
    src_name=strcat(num2str(angle(ag)),'d_',num2str(distance(ds)*100),'cm');
    disp(sprintf('Calculating impulse response for source %s',char(src_name)));
    [H_OUT LEAD_ZEROS] = room_impulse(sources(i,:),recs,walls,wtypes,sph_cent,sph_r,f_samp,c_snd,num_taps,highpass,dsply,err);
    eval(sprintf('save human_hrtf_full_%s *',src_name));
    clear H_OUT LEAD_ZEROS
    end
end
