%
% main_create_stimuli.m
%
% computes one HRTF using spherical/image model
%
% Goal: creates 2-receiver bilateral room impulse responses using a
%       modified version of the room-image method that does not allow for
%       non-integer delays
%

addpath('slamam\modSphericalHRTF');   
addpath('../');   



%% Inline funxtions
% Signal normalization; avoide clipping when saving to WAV file
sig_norm = @(x) x./max(abs(x));


%% Folders & Generic File Names
% the path to where the WAV files will be saved 
path2save  = '.\.stim2Impal\';        

% the path to TIMIT database; this is the source of the DRY stimulus
path2TIMIT = 'D:\DataBases\TIMIT';  
% path2TIMIT = '.\.database\';  

% The new file name to save
% fn_dry = 'TIMIT_dry';


%% The stimuli sampling rate to load\use
Fs = 100e3;     % (Hz)
%Fs = 16e3;      % (Hz)


%% Load Concatenated Stimulus

% Load table TBL_IMF; this table contains IMF fields of the TIMIT sentences 
path_root = load.path_to_data('project1');
load( [path_root, 'Spch_(36)sec_new', filesep, 'tbl_imf.mat'] );

% Save the desired IMFs into 
selected_imf = 2;
fprintf('--> Selected IMF: %d\n', selected_imf);

path2imfs = [path2save, sprintf('IMF%d\\', selected_imf)];
mkdir(path2imfs);

% Save desired IMFs of all sentences
for q = 1:size(tbl_imf,1)
    imf2save = tbl_imf.imf{q}(:, selected_imf);
    
    % Upsample to Fs (if needed)
    fs_load = tbl_imf.fs_short(q);
    [num, den] = rat(Fs/fs_load);
    imf2save = resample(imf2save, num, den);
    
    filename_imf = sprintf('%simf(%d)_%02d.wav', path2imfs, selected_imf, q);
    
    %save(filename_imf, 'imf2save');
    audiowrite(filename_imf, imf2save, Fs);

end
clear imf2save



%% Concatenate all IMFs into one sentence
% % Create a new stimulus?
% new_stimulus = true;
% 
% % % The stimuli sampling rate to load\use
% % Fs = 100e3;     % (Hz)
% % %Fs = 16e3;      % (Hz)
% 
% switch Fs
%     case 100e3
%         path_females = load.path_to_data('timit_females_36sec_fs100khz');
%         %path_males   = load.path_to_data('timit_males_36sec_fs100khz');
% 
%     case 16e3
%         path_females = load.path_to_data('timit_female');
%         %path_males   = load.path_to_data('timit_male');
% 
%     otherwise
%         error('--> ERROR: unrecognized Fs!');
% end
% 
% % Check that the paths are valid
% assert(~isempty(dir(path_females)), '--> <path_females> is invalid!');
% assert(~isempty(dir(path_males)), '--> <path_male> is invalid!');

% Concatenate the signals
[y_cat, wav_file_list, wav_start, wav_end, min_silence_sec, min_silence_smp] = ...
    concatenate_wav(...
        'Fs', Fs,...
        'path_females', path2imfs, ...
        'path_males', '' ...
    );

%assert(0 == seg_time_samples, '--> ERROR: set <seg_time_samples> to be zero!');
% duration_sec = len_yv/Fs;  % (sec) duration of the concatenated signal


%% Save the new aggregated audio file for later reference   
fprintf('--> Saving the new aggregated audio file...\n');
filename_cat = sprintf('%sTIMIT_IMF(%d)_concatenate.wav', path2save, selected_imf);
audiowrite(filename_cat, y_cat, Fs);
    



%% ====== RIR ======
addpath('slamam\modSphericalHRTF');   

% Sampling rate
Fs = 100e3;     % (Hz)

% The room's size
walls = [11 13 3];      % (m)

% REVERBERATION: the walls type (see acoeff.m for more details)
% Project 1:  [1-10*eps, 0.8, 0.2]
% Project 2:  [1-10*eps, 0.75, 0.44, 0.11]
wtypes = 0.11;  % [1-10*eps, 0.75, 0.44, 0.11];   %[1-10*eps, 0.5, 0.2];

% Source(s)
angle = 0;      % [-90:15:90];  	% degrees
dist  = 3.0;    % [0.5 1 1.5 3];    

% SPHERE center 
sph_cent  = [4.7 5.4 1.4];   % (m) sphere is in center of room
recs_dist = 10.29 * 0.01;    % (m) rabbit: 10.29 cm    

% The length of the RIR 
rir_len_s = 3.0;    % (sec) 
num_taps  = ceil(rir_len_s * Fs); % (smp) # of samples to calculate for the RIR

% Calculate the RIR
tbl_rir = RIR(num_taps, Fs, walls, wtypes, angle, dist, sph_cent, recs_dist);

% Add row numbers to the table
tbl_rir = [array2table((1:size(tbl_rir,1))'), tbl_rir];
tbl_rir.Properties.VariableNames{1} = 'row';

% Reorder the columns of the table (move the DRRs columns forward)
tbl_rir = tbl_rir(:,[1:4, 8:9, 5:7, 10:end]);
    
disp(tbl_rir);



%% Convolve with the reverberation RIR
y_rir = ConvRIR(y_cat, tbl_rir.rir);





%% Split & SAVE the WAV files for Impale
% fn_generic = 'spch';
% max_seg_time = 1;   % (sec)


% % Choose the columns for the outer and inner loops in Impale:
% outer_inner_cols = [2, 1];

opts.Fs             = Fs;
opts.fn_generic   	= sprintf('spch_imf(%d)', selected_imf);
opts.rms          	= 'win';  
opts.path2save    	= path2save;
opts.duration_sec 	= 1.0;   % sec          % max segment length to split, in seconds, of each interval (for Impale)
opts.params.name 	= {'Dist', 'Rev'};      % list of all parameters to add to teh file name
opts.params.val     = tbl_rir{:, [2, 1]};   % (dist, wtypes, angle)
% opts.params.name 	= {'dst', 'rev', 'az'}; % list of all parameters to add to teh file name
% opts.params.val     = tbl_rir{:, [3, 2, 4]};   % (dist, wtypes, angle)
opts.params.val(:,1)= ceil(100 * opts.params.val(:,1));    % (cm) Distance in cm
opts.params.val(:,2)= ceil(100 * opts.params.val(:,2));    % (%) Absorption coeff, in percents

%rms_values = split_wav_to_Impale( Yv, Fs, tbl_rir, outer_inner_cols, fn_generic, path2save, max_seg_time );
[rms_values, file_list] = split_wav_to_Impale( y_rir, tbl_rir, opts );




%% Create the WAV information table (the "meta"-information)
% Important -- save the actual samples of the starting and ending of each
% segment!
tbl_tmp = table(wav_start, wav_end);

if exist('tbl_data', 'var') && ~isempty(tbl_data)
    aux.cprintf('cyan', '--> using LOADED <tbl_data>!!\n');
    
    tbl_data.wav_start = tbl_tmp.wav_start;
    tbl_data.wav_end = tbl_tmp.wav_end;
else
    %error('--> Currently doesn''t work; Use main_create_Tmeta.m instead!');
    tbl_data = extract_TIMIT_meta( wav_file_list, path2TIMIT );
    
    % Re-arrange the columns
    tbl_data = [tbl_data(:,1), tbl_tmp, tbl_data(:,2:end)];    
end

clear tbl_tmp



%% Save the meta-data of the stimuli
% dry_wav_fn = [path2save, filesep, '_', fn_dry, '.wav'];
% audiowrite( dry_wav_fn, yv, Fs );
% aux.cprintf('cyan', '--> Dry WAV file was saved (%s)\n', dry_wav_fn);

meta_fn = [path2save, '_spch_metadata.mat'];
save( meta_fn, 'tbl_rir', 'tbl_data', 'Fs', 'file_list', 'rms_values',...
    'min_silence_sec', 'min_silence_smp' );
aux.cprintf('cyan', '--> META-data file was saved (%s)\n', meta_fn);












