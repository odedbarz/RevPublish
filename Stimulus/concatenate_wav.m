function [yv, fullfile_list, wav_start, wav_end, min_silence_sec, min_silence_smp] = ...
    concatenate_wav( varargin )
%
% function [yv, fullfile_list, wav_start, wav_end] = concatenate_wav( varargin )
%
% Description:
% * Adds (aggregates) together all wav files into one big file;
% * Each WAV file is seperated by silence.
%

%% Parse the input
p = inputParser;

addOptional(p, 'Fs', 16e3, @isnumeric);
addOptional(p, 'min_silence_sec', 0.2, @isnumeric);     % 0.2 sec
addOptional(p, 'path_females', load.path_to_data('timit_females_36sec_fs100khz'), @ischar);
addOptional(p, 'path_males', load.path_to_data('timit_males_36sec_fs100khz'), @ischar);
addOptional(p, 'perm_list', 1, @(x) isnumeric(x) || islogical(x));
addOptional(p, 'fignum', [], @isnumeric);

parse(p, varargin{:});

Fs              = p.Results.Fs;
min_silence_sec = p.Results.min_silence_sec;
path_females    = p.Results.path_females;
path_males      = p.Results.path_males;
perm_list       = p.Results.perm_list;
fignum          = p.Results.fignum;



%%
% Minimum size of each slice
min_silence_smp = ceil(Fs*min_silence_sec);
if 0 == min_silence_sec
    fprintf('--> [concatenate_wav.m]: There is NO ZERO PADDING between the TIMIT sentences!');
end

% Initiate the file list to aggregate
files = [];

% The females
path_females = fullfile(path_females, filesep);     % make sure that the last char is a file separeto
files.females = arrayfun(@(S) [path_females, S.name], dir(path_females), 'UniformOutput', 0);
files.females = files.females(3:end);

% The males
path_males = fullfile(path_males, filesep);     % make sure that the last char is a file separetor
files.males   = arrayfun(@(S) [path_males, S.name], dir(path_males), 'UniformOutput', 0);
files.males   = files.males(3:end);

fullfile_list = [files.females; files.males];
n_files        = length(fullfile_list);

% The new audio aggregated file
yv = [];

% Randomize the list
if 1 == perm_list
    Iperm = randperm( n_files );
    %Iperm = [8, 10, 11, 6, 9, 4, 5, 1, 7, 3, 2, 12]; '[concatenate_wav.m]: *** USING FIXED WAV FILES PERMUTATION ***'
    fprintf('--> [concatenate_wav.m]: ### Overriding the permutation of the WAV sequances ###\n');
    fullfile_list = fullfile_list(Iperm);
end

%%
% Saves the starting & ending of each WAV file in the concatenated file.
%seg_time_samples = zeros(Nfiles, 2);
wav_start = zeros(n_files, 1);
wav_end   = zeros(n_files, 1);

seg_prev_start = 0;     % init.
for kk = 1:n_files
    % Full filename (including path)
    fn_kk = fullfile_list{kk};    
    
    [y_load, Fs_load] = audioread( fn_kk );
    assert(Fs == Fs_load,...
        sprintf('--> ERROR: All loaded audio files MUST have Fs=%g Hz sampling rate!', Fs));
    
    % Normalize the audio wav file
    y_load = 0.999*y_load/max(abs(y_load));
    
    % Add the audio file
    % * pick just ONE CHANNEL;
    % * pad with zeros; 
    % * make sure that the minimum distance between two segments is
    %   MIN_SILENCE_SMP;
    y_seg_kk = [y_load(:,1); zeros(min_silence_smp, 1)];
    
    % Add to the overall signal
    yv = [yv; y_seg_kk];  
    
    %seg_time_samples(kk,:) = seg_prev_start + [1, size(y_seg_kk, 1)];
    wav_start(kk) = seg_prev_start + 1;
    wav_end(kk)   = seg_prev_start + size(y_seg_kk, 1);
    
    %seg_prev_start = seg_time_samples(kk,2);
    seg_prev_start = wav_end(kk);
end



%% Pad with zeros at the end
len_yv = length(yv);           % (smp)
duration_sec = len_yv/Fs;    % (sec) signal duration before zero padding

new_duration_sec = ceil(duration_sec);
n_zero_padding = round((new_duration_sec-duration_sec)*Fs);

yv = [yv(:); zeros(n_zero_padding, 1)];
len_yv = length(yv);           % (smp)

% Make sure that the length of the padded signal is OK
% * Note: usually, it's unnecessary to use ceil\round\fix in the assert
assert(  new_duration_sec == len_yv/Fs, '-> ERROR: something is wrong with the zero padding!!' );


fprintf('-> [concatenate_wav.m]:\n');
fprintf('---> Total time: %.2f sec (%.2f min)\n', duration_sec, duration_sec/60);
fprintf('---> Nyv       : %d samples\n', len_yv);

if ~isempty(fignum)
    tt = linspace(0, 1/Fs*size(yv,1), size(yv,1));

    % Plot the new signal
    figure(fignum);
    clf;
    plot(tt, yv);
    xlabel('Time (sec)');
    ylabel('Amplitude');
end















