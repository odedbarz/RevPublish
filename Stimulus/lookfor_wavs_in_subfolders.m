function n_files = lookfor_wavs_in_subfolders(path_corpus, path_subs, save_2_path, n_files)
%
% function n_files = lookfor_wavs_in_subfolders(path_corpus, path_subs, save_2_path, [n_files])
%
% Input:
%   path_corpus: (str) absolute path to the TIMIT corpus;
%   path_subs  : (str) relative path to the next subfolder;
%   save_2_path: (str) a path to save all the WAV files that the function
%                finds in the subfolders.
%
%
% Description:
% Extract TIMIT wav files from the TIMIT corpus. The filenames of the wav files 
% are composed of their path location in the TIMIT corpus, so it's easy to
% locate them and their associated metadata (phones, words, etc).
%
%
% From the TIMIT help function:
%
%     4. CDROM TIMIT Directory and File Structure
%     -- ----------------------------------------
% 
%     The speech and associated data is organized on the CD-ROM according to the
%     following hierarchy:
% 
%     /<CORPUS>/<USAGE>/<DIALECT>/<SEX><SPEAKER_ID>/<SENTENCE_ID>.<FILE_TYPE>
% 
%          where,
% 
%          CORPUS :== timit
%          USAGE :== train | test
%          DIALECT :== dr1 | dr2 | dr3 | dr4 | dr5 | dr6 | dr7 | dr8 
%                      (see Table 1 for dialect code description)
%          SEX :== m | f
%          SPEAKER_ID :== <INITIALS><DIGIT>
% 
%               where, 
%               INITIALS :== speaker initials, 3 letters
%               DIGIT :== number 0-9 to differentiate speakers with identical
%                         initials
% 
%          SENTENCE_ID :== <TEXT_TYPE><SENTENCE_NUMBER>
% 
%               where,
% 
%               TEXT_TYPE :== sa | si | sx
%                             (see Section 2 for sentence text type description)
%               SENTENCE_NUMBER :== 1 ... 2342
% 
%          FILE_TYPE :== wav | txt | wrd | phn
%                        (see Table 5 for file type description)
% 
%     Examples:
%          /timit/train/dr1/fcjf0/sa1.wav
% 
%          (TIMIT corpus, training set, dialect region 1, female speaker, 
%           speaker-ID "cjf0", sentence text "sa1", speech waveform file)
% 
%           /timit/test/df5/mbpm0/sx407.phn
% 
%           (TIMIT corpus, test set, dialect region 5, male speaker, speaker-ID
%            "bpm0", sentence text "sx407", phonetic transcription file)
% 
%

if 4 > nargin
    n_files = [Inf, Inf];  % (1x2) (# Females, # Males) load all available files
elseif 1 == numel(n_files)
    n_files = n_files*[1, 1];
end

new_dir = dir([path_corpus, path_subs]);

for k = 1:length(new_dir)
    if strcmp(new_dir(k).name, '.') || strcmp(new_dir(k).name, '..')
        continue;
    end

    % Recirsive part
    if 1 == new_dir(k).isdir
        next_level_path_2_look = fullfile(path_subs, new_dir(k).name);
        n_files = lookfor_wavs_in_subfolders(path_corpus, next_level_path_2_look, save_2_path, n_files);
    end
    
    % Get only the WAV files
    fn = new_dir(k).name;
    [~, fileparts_name, fileparts_ext] = fileparts(fn);
    if ~strcmpi('.wav', fileparts_ext)
        continue;
    end
    
    % Save the WAV file
    fn = fullfile(path_subs, fileparts_name);
    fn_full = fullfile(path_corpus, [fn, fileparts_ext]);
            
    
    % Get the speaker's gender
    file_props = regexp(fn, filesep, 'split');
    gender = file_props{4}(1);  % 'F' for female or 'M' for male
    if strcmpi('F', gender)
        gender_subfolder = 'Females';        
    elseif strcmpi('M', gender)
        gender_subfolder = 'Males';        
    else
        error('--> ERROR in [lookfor_wavs_in_subfolders.m]: gender MUST be either ''F'' for female or ''M'' for male');
    end

    % This code was taken from the TIMITreagGUI.m GUI
    %%
    fp = fopen(fn_full,'r');
    if fp == -1 ; error('Error on opening file'); return; end;

    fseek(fp,0,'bof');
    nist=fscanf(fp,'%s',1);
    if strcmp(nist(1:4),'NIST') ==0
        %error('Error reading voice file header, not TIMIT ...');
        warning('\n--> I can''t read <%s>\n--> Skipping this file!', fn_full)
        continue;
    end

    % number of bytes for file header
    nbyte_header=fscanf(fp,'%d',1);

    % skip database_id, database_version, utterance_id,
    % channel_count,sample_count
    for i=1:5; X=fscanf(fp,'%s %s %s',3); end

    % obtain sample rate
    X=fscanf(fp,'%s %s',2);  SampleRate =fscanf(fp,'%d',1);

    % skip sample_max and sample_min
    for i=1:2; X=fscanf(fp,'%s %s %s',3); end;

    % obtain bytes per sample
    X=fscanf(fp,'%s %s',2); bps=fscanf(fp,'%d',1);
    if bps == 2;   ftype='short'; end

    % skip the header
    fseek(fp,nbyte_header,'bof');

    wave = fread(fp,inf,ftype);
    fclose(fp);

    maxmag = max(abs(wave));
    wave = 0.95*wave/maxmag;    
    

    %% Set the filename & save the WAV     
    filename = sprintf('%s_', string(file_props) );
    filename = filename(1:(end-1));         % remove the last '_'
    
    
    if strcmpi('F', gender)
        if 0 >= n_files(1)
            return;     % move on to another directory (the other sex)
        end
        n_files(1) = n_files(1) - 1;
        
    elseif strcmpi('M', gender)
        if 0 >= n_files(2)
            return;     % move on to another directory (the other sex)
        end
        n_files(2) = n_files(2) - 1;
        
    else
        error('--> ERROR in [lookfor_wavs_in_subfolders.m]: gender MUST be either ''F'' for female or ''M'' for male');
    end

    
    % Save the wav file
    path_k = fullfile(save_2_path, gender_subfolder);
    if ~isdir(path_k)
        mkdir(path_k);
    end
    
    fn_2_save = fullfile(path_k, [filename, fileparts_ext]);    
    audiowrite(fn_2_save, wave, SampleRate);


end





