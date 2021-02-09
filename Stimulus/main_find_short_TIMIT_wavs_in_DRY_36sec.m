%
% main_find_short_TIMIT_wavs_in_DRY_36sec.m
%
% Description:
% A simple script to find the starting and ending indices of the TIMIT
% sentences in the concatenated dry sentence. Creates a table
% (tbl_metadata) that holds all of this information.
% 

clc

STIM_TYPE_STR = '36';   % the stimulus length

% Load the Tmeta table
fn.path = sprintf('%s%s%s',...
    'C:\Users\barzelo\Google Drive\codeOnCloud\Reverberation\Stimulus\.database\Project1\Spch_(', ...
    STIM_TYPE_STR, ...
    ')sec\');
fn.file  = sprintf('spch_%s_metadata.mat', STIM_TYPE_STR);
load( fullfile(fn.path, fn.file) );

% Get paths to files
files.timit_short = {'timit_females_36sec', 'timit_males_36sec'};  
files.timit_full = {'wav_female_36sec_full_length', 'wav_male_36sec_full_length'};


%% Create the WAV table
tbl_wav = Tmeta(:,1);       % the 'fn' column
n_rows  = size(tbl_wav,1);   % == 12
tbl_wav = [tbl_wav, array2table( zeros(n_rows,2), 'VariableNames', {'t0', 't1'} )];
tbl_wav = [tbl_wav, array2table( zeros(n_rows,2), 'VariableNames', {'t0_timit', 't1_timit'} )];
tbl_wav = [tbl_wav, array2table( cell(n_rows,1), 'VariableNames', {'wav_short'} )];
tbl_wav = [tbl_wav, array2table( cell(n_rows,1), 'VariableNames', {'wav_full'} )];
tbl_wav = [tbl_wav, array2table( zeros(n_rows,1), 'VariableNames', {'fs'} )];


%% Load the concatenated sentence
[path_stim, ~] = load.path_to_data('wav_spch_36sec');    
[wav_stim, fs_stim] = audioread([path_stim, 'BRIR_L_-Dist15_-Rev10.wav']);

% Resample to 16kHz (if needed)
fs = 16e3;  % (Hz) sampling rate
[num, den] = rat(fs/fs_stim);
wav_stim = resample(wav_stim, num, den);



%% load all WAVs
fprintf('\n');
fprintf('--> loading all WAV files...\n');

files.fullpath = {};
files.fn = {};

for m = 1:n_rows  
    % Get the m'th TIMIT wav file
    file_m = Tmeta.fn{m};
    idx_m = m;

    if contains(file_m, '_F')       % is a male speaker?
        timit_full = files.timit_full{1};
        timit_short = files.timit_short{1};
    elseif contains(file_m, '_M')   % is a female speaker? 
        timit_full = files.timit_full{2};
        timit_short = files.timit_short{2};
    else
        error('--> Something is dead wrong!');
    end
    [path_full,  ~] = load.path_to_data( timit_full );
    [path_short, ~] = load.path_to_data( timit_short );


    %% Load the SHORT TIMIT signal
    [wav_short, fs_short] = audioread([path_short, file_m, '.wav']);

    % Resample to 16kHz (if needed)
    [num, den] = rat(fs_short/fs);
    wav_short = resample(wav_short, num, den);

    tbl_wav.wav_short{idx_m} = wav_short;


    %% Load the FULL TIMIT signal
    [wav_full, fs_full] = audioread([path_full, file_m, '.wav']);

    % Resample to 16kHz (if needed)
    [num, den] = rat(fs_full/fs);
    wav_full = resample(wav_full, num, den);

    tbl_wav.wav_full{idx_m} = wav_full;


    %% Set the sampling frequency
    tbl_wav.fs(idx_m) = fs;

end
fprintf('--> All files were loaded!\n');



%% Find the starting & ending points 
fprintf('\n');
fprintf('--> Find the starting & ending points...\n');

y = flipud(wav_stim);

for k = 1:size(tbl_wav,1)
    % Find the starting point of the short WAV in the full WAV    
    [C, lags] = xcorr( tbl_wav.wav_full{k}, tbl_wav.wav_short{k} ); % (Option #1) cross-correlation
    [~, idx_max] = max(C);       
    tbl_wav.t0_timit(k) = lags(idx_max) +1;   % (smp) starting point of the short TIMIT sentence     
    tbl_wav.t1_timit(k) = tbl_wav.t0_timit(k) + length(tbl_wav.wav_short{k}) -1; 
    
    
    % Find the starting point of the short WAV in the WHOLE stimulus (i.e., the concatenated WAV file)
    x = tbl_wav.wav_short{k};   
    C_ = fftfilt(x, y);   % (Option #2) cross-correlation
    [~, idx_max] = max(C_);    
    tbl_wav.t0(k) = length(y) - idx_max +1;

    % (smp) ending point of the short TIMIT sentence
    tbl_wav.t1(k) = tbl_wav.t0(k) + length( tbl_wav.wav_short{k} ) -1;
    
    % Make sure that these segments have the same length!
    assert(tbl_wav.t1_timit(k) - tbl_wav.t0_timit(k) ==  tbl_wav.t1(k) - tbl_wav.t0(k));
end



%% Test
figure(999);
clf;

idx2show = 5;

ax = subplot(2,1,1);
x = ( 1:length(wav_stim) )';
x_short = x(tbl_wav.t0(idx2show):tbl_wav.t1(idx2show));
plot( x, wav_stim );
hold on
plot(x_short , tbl_wav.wav_short{idx2show} );
hold off
title( aux.mName2latex(tbl_wav.fn{idx2show}) );
legend('Full Stimulus', sprintf('TIMIT short wav'));
axis tight
ylabel('Full');

ax(2) = subplot(2,1,2);
x_ = ( 1:length(tbl_wav.wav_full{idx2show}) )';
x_short_ = x_(tbl_wav.t0_timit(idx2show):tbl_wav.t1_timit(idx2show));
plot( x_, tbl_wav.wav_full{idx2show} );
hold on
plot( x_short_, tbl_wav.wav_short{idx2show} );
hold off
axis tight
disp(tbl_wav);
ylabel('TIMIT WAV');
xlabel('Samples');


%% 
info = '* t0 & t1 gives the starting and ending times in which each TIMIT WAV file starts in the main stimulus.';
info = [info, '\n* t0_timit & t1_timit gives the starting and ending times for each short WAV in the full WAV. Use these times to get the phonemes & words of the original TIMIT WAV files.\n'];
fprintf(info);
fprintf('\n');
tbl_metadata = join(tbl_wav, Tmeta);


%% Save everything into one table
%{
    '### SAVE the metadata table! ###'
    file_2_save = sprintf('../.data/metadata_(%s)_wav_(%s).mat', STIM_TYPE_STR, date);
    fprintf('-> Saving the METADATA into file <%s>\n', file_2_save);
    save(file_2_save, 'tbl_metadata');
%}






