%
% main_white_noise_interferer_scipt.m.m
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

% Inline funxtions
% Signal normalization; avoide clipping when saving to WAV file
sig_norm = @(x) x./max(abs(x));

% The path to where the WAV files will be saved 
path2save  = '.\.stim2Impal\';        

Fs = 100e3;                 % (Hz)
duration_sec = 3.0;     % (sec) stimulus length of all segments



%% Create TARGET Stimuli
% path_to_wav_files{1} = load.path_to_data('project2_timit_females');
% path_to_wav_files{2} = load.path_to_data('project2_timit_males');
% itr = main_generate_interferer_stimuli(path_to_wav_files, Fs, duration_sec, []);
itr.Fs = Fs;
itr.duration_init_sec = 2 * duration_sec;     % (sec) stimulus length of all segments
itr.duration_init_smp = fix(itr.duration_init_sec * itr.Fs);     % (sec) stimulus length of all segments
itr.ywn_all = randn(itr.duration_init_smp, 1);



%% Trir: INTERFERER's RIR
new_rir = 1;

if 1 == new_rir
    main_RIR_script;
else
    % Load TRIR table
    load('.\.database\tbl_RIR_Dist(3m).mat');     % --> Trir
end



%% Convolution --- Interferer
% Cyclic convolution
% Ywn_all = fft(itr.ywn_all);
% yrir_all = cellfun(@(H) ifft(fft(H).*Ywn_all, 'symmetric'), Trir.rir, 'UniformOutput', false);
yrir_all = cellfun(@(B) fftfilt(B, itr.ywn_all), Trir.rir, 'UniformOutput', false);

% Remove the build-up of the reverberation (after the stimulus' onset)
duration_smp = fix(duration_sec * itr.Fs);  % sampples
yrir = cellfun(@(Y) Y(end-[(duration_smp-1):-1:0],:), yrir_all, 'UniformOutput', false);
% yrir = cell(1, length(yrir_all));
% for k = 1:length(yrir_all)
%     yrir{k} =  yrir_all{k}(end-[(duration_smp-1):-1:0], :);
% end


% The stimulus' table
itr.Tstim = [Trir(:,1:6), table(yrir)];



%% The interaural coherence
coherence = interaural_coherence(itr.Tstim.yrir);
itr.Tstim = [itr.Tstim, table(coherence)];



%% SAVE everything --- Interferer wav files ---
opts.Fs                  = Fs;     % (Hz)
opts.path2save           = '.stim2Impal/';
opts.duration_sec        = duration_sec;   % sec

itr.opts                 = opts;
itr.opts.fn_generic      = 'WhiteNoise_itr'; 
itr.opts.rms             = 'simple';

itr.opts.params.name     = {'dst', 'rev', 'az'};
itr.opts.params.val      = itr.Tstim{:, [3, 2, 4]};
itr.opts.params.val(:,1) = ceil(100 * itr.opts.params.val(:,1));    % Distance in cm
itr.opts.params.val(:,2) = ceil(100 * itr.opts.params.val(:,2));    % Absorption coeff in %


itr.rms = split_wav_to_Impale( itr.Tstim.yrir, Trir, itr.opts );

%% Save metadata into XLS file 
writetable( Trir(:,{'row', 'wtypes', 'dist', 'angle', 'drrL', 'drrR'}), ...
    [itr.opts.path2save, '!metadata.xls'], 'Sheet', 'RIR');

col_not_2_include = find(strcmpi(itr.Tstim.Properties.VariableNames, 'yrir'));
writetable( [itr.Tstim(:,1:end~=col_not_2_include), array2table(itr.rms, 'VariableNames', {'rms_L', 'rms_R'})], ...
    [itr.opts.path2save, '!metadata.xls'], 'Sheet', 'ITR');


%% Save metadata into MAT file 
save([itr.opts.path2save, '!itr_data.mat'], 'itr', 'Trir');





