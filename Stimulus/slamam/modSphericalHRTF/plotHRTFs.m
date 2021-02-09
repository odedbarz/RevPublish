%%MAKE SURE TO CHANGE SAMPLING RATE!!%
addpath c:\delgutte\devores\catHRTF\HRTF_FULL\CleanHRTFs\Normalized

dis = [50 100 150 ];
ang = [-90 -75 -60 -45 -30 -15 -4 15 30 45 60 75 90];

SAMPF=50000

for d_ind=3:3
    figure;
    for a_ind = 1:length(ang)
        eval(sprintf('load anech_full_%dd_%dcm.mat', ang(a_ind), dis(d_ind)));
        subplot(4,4,a_ind);
        plot((1:length(wimpl))/SAMPF,wimpl); hold on; plot((1:length(wimpr))/SAMPF,wimpr,'r');
        xlabel('time (s)');
        title(sprintf('Anechoic; %d deg.azimuth, %d cm distance',ang(a_ind),dis(d_ind)));
    end
    legend('left', 'right');
end


for d_ind=1:3
    figure;
    for a_ind = 1:length(ang)
        eval(sprintf('load reverb_full_%dd_%dcm.mat', ang(a_ind), dis(d_ind)));
        subplot(4,4,a_ind);
        plot((1:length(fimpl))/SAMPF,fimpl); hold on; plot((1:length(fimpr))/SAMPF,fimpr,'r');
        xlabel('time (ms)');
        title(sprintf('Reverberant; %d deg. azimuth; %d cm distance',ang(a_ind),dis(d_ind)));
    end
    legend('left', 'right');
end

