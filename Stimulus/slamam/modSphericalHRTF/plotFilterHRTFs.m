%%MAKE SURE TO CHANGE SAMPLING RATE!!%

dis = [50 100 150 300];
ang = [-90 -75 -60 -45 -30 -15 0 15 30 45 60 75 90];

SAMPF=50000

for d_ind=1:3
    figure;
    for a_ind = 1:length(ang)
        eval(sprintf('load anech-azitd%d_%dcm_left.mat', ang(a_ind), dis(d_ind)));
        subplot(4,4,a_ind);
        plot((1:length(filt.b))/SAMPF,filt.b); hold on; 
        eval(sprintf('load anech-azitd%d_%dcm_right.mat', ang(a_ind), dis(d_ind)));
        plot((1:length(filt.b))/SAMPF,filt.b,'r');  
        xlabel('time (s)');
        title(sprintf('Anechoic; %d deg.azimuth, %d cm distance',ang(a_ind),dis(d_ind)));
    end
    legend('left', 'right');
end


for d_ind=1:3
    figure;
    for a_ind = 1:length(ang)
        eval(sprintf('load reverb-azitd%d_%dcm_left.mat', ang(a_ind), dis(d_ind)));
        subplot(4,4,a_ind);
        plot((1:length(filt.b))/SAMPF,filt.b); hold on; 
        eval(sprintf('load reverb-azitd%d_%dcm_right.mat', ang(a_ind), dis(d_ind)));
        plot((1:length(filt.b))/SAMPF,filt.b,'r');  
        xlabel('time (s)');
        title(sprintf('Reverberant; %d deg. azimuth; %d cm distance',ang(a_ind),dis(d_ind)));
    end
    legend('left', 'right');
end

