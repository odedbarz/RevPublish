% 
%   main_RIR_script.m
%
%      [  1 | 0.5 |  0.2 |   0.1]
% DRR: [Inf | 3.6 | -2.5 | -11.6]
%

addpath('../');
addpath('./slamam/modSphericalHRTF');   

% Sampling rate
Fs = 100e3;     % (Hz)

% The room's size
walls = [11 13 3];      % (m)

% REVERBERATION: the walls type (see acoeff.m for more details)
% Project 1:  [1-10*eps, 0.5, 0.2]
% Project 2:  [1-10*eps, 0.75, 0.44, 0.11]
wtypes = [1-10*eps, 0.8, 0.2];     %[1-10*eps, 0.75, 0.44, 0.11];   

% Source(s)
angle = 0;      %-90:15:90;  % degrees
dist  = [1.5, 3.0];                 % (m) distance between source and target; def: [0.5 1 1.5 3];    

% SPHERE center 
sph_cent  = [4.7 5.4 1.4];   % (m) sphere is in center of room
recs_dist = 10.29 * 0.01;    % (m) rabbit: 10.29 cm    

% The length of the RIR 
rir_len_s = 3.0;    % (sec) 
num_taps  = ceil(rir_len_s * Fs); % (smp) # of samples to calculate for the RIR

% Calculate the RIR
Trir = RIR(num_taps, Fs, walls, wtypes, angle, dist, sph_cent, recs_dist);

% Add row numbers to the table
Trir = [array2table((1:size(Trir,1))'), Trir];
Trir.Properties.VariableNames{1} = 'row';

% Reorder the columns of the table (move the DRRs columns forward)
Trir = Trir(:,[1:4, 8:9, 5:7, 10:end]);




% '*** SAVE !! ***'
% save('./.stim2Impal/tbl_RIR_Dist(3m).mat', 'Trir')

% round(10*mean([Trir.drrL, Trir.drrR],2))/10

%% Plot
fignum = 1
if ~isempty(fignum)
    figure(fignum);
    clf;
    
    idx2plot = 2;
    
    t = linspace(0, num_taps/Fs, num_taps)';
    warning off
    plot(t, 20*log10(Trir.rir{idx2plot}));
    warning on
    ylabel('RIR $(20log_{10}(\cdot))$');
    xlabel('Time (sec)');
    
    title_1 = sprintf('Room Impulse Response');
    title_2 = sprintf('DRR$_L$: %.1f dB, DRR$_R$: %.2f dB; Dist: %g m; Revb: %g',...
        Trir.drrL(idx2plot), Trir.drrR(idx2plot), ...
        Trir.dist(idx2plot), Trir.wtypes(idx2plot) );
    title(sprintf('\\begin{tabular}{c} %s \\\\ %s \\end{tabular}', title_1, title_2), 'Interpreter', 'latex');
    
    legend('Left', 'Righ');
    ylim([-100, 0]);
    xlim([0, t(end)]);
    hold on
    line_h = line(xlim, -60*[1 1]);
    line_h.Color = 'k';
    line_h.LineStyle = '--';
    line_h.LineWidth = 2;
    text_h = text(0.9*max(xlim), -60, '$T_{60}$');
    text_h.Interpreter = 'latex'; 
    text_h.FontSize = 32; 
    text_h.VerticalAlignment = 'bottom'; 
    hold off
end



