%
% main_target_scipt.m
%
% Project #2 (Bertrand & Joe).
%


clc
close all
clear all

addpath('slamam\modSphericalHRTF');   
addpath('..');   

set(0,'defaultAxesFontSize',18)
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');
set(0,'defaultTextInterpreter','latex'); %trying to set the default

% The path to where the WAV files will be saved 
path2save  = '.\.stim2Impal\';        

Fs = 100e3;         % (Hz)
duration_sec = 3;    % (sec) stimulus length of all segments


% gendre_str = 'F';
gendre_str = 'M';



%% Create TARGET Stimuli
trg = main_generate_target_stimuli(Fs, gendre_str, duration_sec, []);



%% TARGET's RIR 
new_rir = 0;

if 1 == new_rir
    main_RIR_script;
else
    % Load TRIR table
    load('.\.database\tbl_RIR_Dist(3m).mat');     % --> Trir
end


%% Convolution --- Target
yrir = cellfun(@(H) fftfilt(H, trg.yspch), Trir.rir, 'UniformOutput', false);
trg.Tstim = [Trir(:,1:4), table(yrir)];



%% The interaural coherence
coherence = interaural_coherence(trg.Tstim.yrir);
trg.Tstim = [trg.Tstim, table(coherence)];



%% SAVE everything --- Write the Target wav files ---
opts.Fs                  = Fs;     % (Hz)
opts.path2save           = '.stim2Impal/';
opts.duration_sec        = duration_sec;   % sec

trg.opts                 = opts;
trg.opts.fn_generic      = sprintf('SA1_%s_fc(%g)Hz_trg', gendre_str, trg.foct); 
trg.opts.rms             = 'win';

trg.opts.params.name     = {'dst', 'rev', 'az'};
trg.opts.params.val      = trg.Tstim{:, [3, 2, 4]};
trg.opts.params.val(:,1) = ceil(100 * trg.opts.params.val(:,1));    % Distance in cm
trg.opts.params.val(:,2) = ceil(100 * trg.opts.params.val(:,2));    % Absorption coeff in %

trg.rms = split_wav_to_Impale( trg.Tstim.yrir, Trir, trg.opts );


%% Save metadata into XLS file 
% %{
writetable( Trir(:,{'row', 'wtypes', 'dist', 'angle', 'drrL', 'drrR'}), ...
    [trg.opts.path2save, '!metadata.xls'], 'Sheet', 'RIR');

writetable( [trg.Tstim(:,1:end~=5), array2table(trg.rms, 'VariableNames', {'rms_L', 'rms_R'})], ...
    [trg.opts.path2save, '!metadata.xls'], 'Sheet', ['TRG ', trg.gendre_str]);



%% Save metadata into MAT file 
save([trg.opts.path2save, '!trg_data_', trg.gendre_str, '.mat'], 'trg', 'Trir');
%}



%% Plot
rev_list = round(100*unique(Trir.wtypes));
disp(rev_list);

fn2plot = [trg.opts.path2save, trg.opts.fn_generic, '_-dst300_-rev%d_-az0L.wav'];
y11  = audioread(sprintf(fn2plot, 11));
y44  = audioread(sprintf(fn2plot, 44));
y75  = audioread(sprintf(fn2plot, 75));
y100 = audioread(sprintf(fn2plot,100));

figure(11);
clf;
t = linspace(0, length(y11)/Fs, length(y11))';
plot(t, [y11, y44, y75 y100] + 2.0*ones(size(y11,1),4)*diag(1:4)-1.0)
xlabel('Time (sec)');
ylabel('Absorption (parameter REV)');
set(gca, 'YTick', (1:2:8), 'YTickLabel', {'11', '44', '75', '100'});

title_str = '%s speaker (azimuth: $0^0$)';
if strcmpi('F', trg.gendre_str)
    title_str = sprintf(title_str, 'Female');
else
    title_str = sprintf(title_str, 'Male');
end

title(title_str);












