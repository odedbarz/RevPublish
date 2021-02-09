%
% main_create_stimuli_without_revb.m
%
%
% Note: use extract_wav_from_TIMIT.m to get the wav files.
%
%

addpath(['..', filesep]);

n_batch = 6;    % # of desired batches
n_files = 6;    % # of files for M\F

% The stimuli sampling rate to load\use
% Fs = 100e3;     % (Hz)
Fs = 16e3;      % (Hz)

% Save for Impale
Fs_2_Impale = 100e3;
[resmp.N, resmp.D] = rat(Fs_2_Impale/Fs);

% 
min_batch_sec = 32;   % (sec)
% min_batch_smp = min_batch_sec*Fs;

% The max file will be:
%   ceil(max_batch_sec + (2*n_files)*0.2) 
%   * ceil() because Impale needs an integer token (for a repetition of 1 sec);
%   * *2 for females & males;
%   * 0.2 is the (default) padding between WAV files.
max_batch_sec = 33;   % (sec)
% max_batch_smp = max_batch_sec*Fs;
fprintf('--> Estimated final WAV files length is: %g\n', ceil(max_batch_sec + (2*n_files)*0.2));


%% Prepare the RIR 
use_current_rir = exist('Trir', 'var') && 1 == size(Trir,1) && Trir.wtypes == 1-10*eps;
if ~use_current_rir
    addpath('slamam\modSphericalHRTF');   

    % Sampling rate
    %Fs = 100e3;     % (Hz)

    % The room's size
    walls = [11 13 3];      % (m)

    % REVERBERATION: the walls type (see acoeff.m for more details)
    % Project 1:  [1-10*eps, 0.8, 0.2]
    % Project 2:  [1-10*eps, 0.75, 0.44, 0.11]
    wtypes = 1-10*eps;  % [1-10*eps, 0.75, 0.44, 0.11];   %[1-10*eps, 0.5, 0.2];

    % Source(s)
    angle = 0;      % [-90:15:90];  	% degrees
    dist  = 1.5;    % [0.5 1 1.5 3];    

    % SPHERE center 
    sph_cent  = [4.7 5.4 1.4];   % (m) sphere is in center of room
    recs_dist = 10.29 * 0.01;    % (m) rabbit: 10.29 cm    

    % The length of the RIR response vector, in units of time and samples 
    rir_len_s = 3.0;    % (sec) 
    num_taps  = ceil(rir_len_s * Fs); % (smp) # of samples to calculate for the RIR

    % Calculate the RIR
    Trir = RIR(num_taps, Fs, walls, wtypes, angle, dist, sph_cent, recs_dist);

    % Add row numbers to the table
    Trir = [array2table((1:size(Trir,1))'), Trir];
    Trir.Properties.VariableNames{1} = 'row';

    % Reorder the columns of the table (move the DRRs columns forward)
    Trir = Trir(:,[1:4, 8:9, 5:7, 10:end]);
end
disp(Trir);




%% Folders & Generic File Names
% path2save  = '.\.stim2Impal\';        
path_root_2_save  = 'D:\DataBases\TIMIT\wav_files\stim_2_Impale\';

% the path to TIMIT database; this is the source of the DRY stimulus
path_2_TIMIT = 'D:\DataBases\TIMIT';  
% path2TIMIT = '.\.database\';  

% The new file name to save
path2files = 'D:\DataBases\TIMIT\wav_files\';
path_all_females = fullfile(path2files, 'Females');
path_all_males = fullfile(path2files, 'Males');

% Check that the paths are valid
assert(~isempty(dir(path_all_females)), '--> <path_all_females> is invalid!');
assert(~isempty(dir(path_all_males)), '--> <path_all_males> is invalid!');

files_all_females = dir(path_all_females);
files_all_females = files_all_females(3:end);   % remove the '.' & '..' entries
[~, I] = sort([files_all_females.bytes]);
files_all_females = files_all_females(I);

files_all_males   = dir(path_all_males);
files_all_males = files_all_males(3:end);       % remove the '.' & '..' entries
[~, I] = sort([files_all_males.bytes]);
files_all_males = files_all_males(I);



%% ===== Get Batch Stimuli =====
% Create the overall concatenated stimulus (also the dry stimulus)
% Get more batches and then select those that are of the same approximate length
Y_list        = cell(1, n_batch);
Y_len         = zeros(1, n_batch);
wav_file_list = cell(1, n_batch);

idx_reservoir.F = 1:length(files_all_females);  % ...*2 for both sexes
idx_reservoir.M = 1:length(files_all_males);

clear fn_idx

% Create 
for k = 1:n_batch
    path_2_save_batch = fullfile(path_root_2_save, ['batch', num2str(k)], filesep);
    if exist(path_2_save_batch, 'dir')
        % Delete old batch folder, if it exists
        rmdir(path_2_save_batch, 's');
    end
    
    % Get N_FILES from the main folder & save them in a new directory
    path_females_k = fullfile(path_2_save_batch, ['F', num2str(k)]);
    if ~exist(path_females_k, 'dir')
        mkdir(path_females_k);
    else
        error('--> ERROR: the folder <%s> already exists! You need to delete it first!', path_females_k);
    end 
    
    path_males_k = fullfile(path_2_save_batch, ['M', num2str(k)]);
    if ~exist(path_males_k, 'dir')
        mkdir(path_males_k);
    else
        error('--> ERROR: the folder <%s> already exists! You need to delete it first!', path_males_k);        
    end
    
    % FEMALE: Select random files from the list; 
    % NOTE: all of this sorting is to get a (more or less) same length
    %       concatenated signals
    OK = true;
    counter = 1;
    fprintf('--> Looking for the right files...\n');    
    fprintf('---> ');    
    while OK
        fprintf('%d ', counter);
        counter = counter + 1;
        
        dummy = randperm(length(idx_reservoir.F), 3* n_files);   
        dummy = sort(dummy);
        dummy = dummy(1:3:end);
        fn_idx.F = idx_reservoir.F(dummy);        
        
        % MALE: Select random files from the list; 
        % NOTE: all of this sorting is to get a (more or less) same length
        %       concatenated signals    
        dummy = randperm(length(idx_reservoir.M), 3* n_files);
        dummy = sort(dummy);
        dummy = dummy(1:3:end);    
        fn_idx.M = idx_reservoir.M(dummy);
        
        batch_duration = 0;     % (sec)
        for q = 1:n_files
            dir_female_q = files_all_females(fn_idx.F(q));
            audio_data = audioinfo( fullfile(dir_female_q.folder, dir_female_q.name) );
            batch_duration = batch_duration + audio_data.Duration;
 
            dir_male_q = files_all_males(fn_idx.M(q));
            audio_data = audioinfo( fullfile(dir_male_q.folder, dir_male_q.name) );
            batch_duration = batch_duration + audio_data.Duration;
        end
        
        if (batch_duration > min_batch_sec) && (batch_duration < max_batch_sec)
            OK = false;
        elseif 50 < counter     % max counter
            OK = false;
            error('\n--> Could not find suitable WAV files\n');
        end
        
    end
    fprintf('\n--> Done!\n\n');    
    
    % Remove selected files (don't repeat in other batches)
    idx_reservoir.F = setdiff(idx_reservoir.F, fn_idx.F);
    idx_reservoir.M = setdiff(idx_reservoir.M, fn_idx.M);
    
    
    % Copy the files
    % ====================================    
    for q = 1:n_files
        %fn_idx = q + (k-1)*n_files;
        
        % Source F\M
        src.F = fullfile( files_all_females(fn_idx.F(q)).folder, ...
            files_all_females(fn_idx.F(q)).name);
        
        src.M = fullfile( files_all_males(fn_idx.M(q)).folder, ...
            files_all_males(fn_idx.M(q)).name);
        
        copyfile( src.F, path_females_k );
        copyfile( src.M, path_males_k );
    end
    
    % Concatenate the signals
    [Y_list{k}, wav_file_list{k}, wav_start, wav_end, min_silence_sec, min_silence_smp] = ...
        concatenate_wav(...
            'Fs', Fs,...
            'path_females', path_females_k, ...
            'path_males', path_males_k, ...
            'perm_list', true ...
        );
    
    Y_len(k) = length(Y_list{k});
    fprintf('\n');
end

max_Y_length = max(Y_len);
max_Y_sec = ceil(1/Fs * max_Y_length);

% [~, I] = sort(Y_len, 'descend');
% Y_list = Y_list(I);     


% Make sure that the ALL_SEGMENT_LENGTH is longer than teh longest
% concatenated signal
fprintf('--> The lengths of all the batches (in sec):\n');
disp( Y_len/Fs );
max_batch_smp = ceil(max_batch_sec + min_silence_sec*(2*n_files))*Fs;
assert( max_Y_sec <= max_batch_smp, '--> ERROR: You need to choose a longer total segment length!');


%% ==== Save the Batch Stimuli ====
% * All stimuli need to have the same length for Impale to work properly!
% * Pad them with zeros.
for k = 1:n_batch
    assert( max_batch_smp >= length(Y_list{k}), '--> This signal is too long!');

    % Pad them with zeros to the maximum length
    n_pad_zeros = max_batch_smp - length(Y_list{k});           
    Yk = [Y_list{k}; zeros(n_pad_zeros, 1)];
   
    % Save the new aggregated audio file for later reference   
    aux.cprintf('cyan', '--> Saving the new aggregated audio file...\n');
    audiowrite([path_2_save_batch, 'spch_cat_batch', num2str(k),'.wav'], Yk, Fs);
    
    
    % Save the metadata of each batch
    % ====================================    
    tbl_data = extract_TIMIT_meta( wav_file_list{k}, path_2_TIMIT );

    % Important -- save the actual samples of the starting and ending of 
    % each segment!
    tbl_tmp = table(wav_start, wav_end);

    % Re-arrange the columns
    tbl_data = [tbl_data(:,1), tbl_tmp, tbl_data(:,2:end)];

    
    
    % Split & Save the stimuli for Impale
    % ====================================
    yRIR = ConvRIR(Yk, Trir.rir);
    
    % For Impale, all stimuli must have a 100kHz sampling rate 
    yk_100k = cellfun(@(X) resample(X, resmp.N, resmp.D), yRIR, 'UniformOutput', false);
    
    opts.Fs             = Fs_2_Impale;
    opts.fn_generic   	= 'spch';
    opts.rms          	= 'win';  
    opts.path2save    	= fullfile(path_root_2_save, ['batch', num2str(k)], filesep);
    opts.duration_sec 	= 1.0;   % sec          % max segment length to split, in seconds, of each interval (for Impale)
    opts.params.name 	= {'batch'};      % list of all parameters to add to teh file name
    opts.params.val     = k;
    assert(0 ~= exist(opts.path2save, 'dir'), '--> ERROR: the folder <%s> SHOULD already exists!', opts.path2save);
    
    %rms_values = split_wav_to_Impale( Yv, Fs, Trir, outer_inner_cols, fn_generic, path2save, max_seg_time );
    [rms_values, file_list] = split_wav_to_Impale( yk_100k, opts );

    % Save the meta file
    meta_fn = fullfile(opts.path2save, ['!_', opts.fn_generic, '_batch', num2str(k), '.mat']);
    save( meta_fn, 'Trir', 'tbl_data', 'Fs', 'file_list', 'rms_values',...
        'min_silence_sec', 'min_silence_smp' );
    aux.cprintf('cyan', '--> META-data file was saved (%s)\n', meta_fn);
    
    fprintf('\n');
    
    
end







