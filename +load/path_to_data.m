function [p, dir_p] = path_to_data(data_type, verbose)
%
%   function p = path_to_data(data_type, verbose)
%
% Input:
%   data_type (str)
%   * 'meas': (default) points to the database direcotry that holds all
%             measurements.
%   * 'wav' : points to the WAV files of the stimulus (reconstruction project; aim 2).
% 
% Description:
%   Returns the full (absolute) path to a desired directory. These directory
%   usually holds measurements, data, WAV files, etc.
%


if 1 > nargin || isempty(data_type)
    data_type = 'meas';
end

if 2 > nargin
    verbose = [];
end

% Main (root) directory
[rootpath, ~, ~] = fileparts( mfilename('fullpath') );

% path_to_data is now a part of the +load library, so get one level "up"
% (i.e., "outside" +load)
[rootpath, ~, ~] = fileparts( rootpath );



switch lower(data_type)        
    case {'data', '_data'}     % case 'meas'/'impale_data'
        p = fullfile(rootpath, '_data', filesep);
        
                
    case 'stimulus'
        path_data = load.path_to_data('_data');
        p = fullfile(path_data, 'Stimulus', filesep);
        
    case 'analysis'
        path_data = load.path_to_data('_data');
        p = fullfile(path_data, 'Analysis', filesep);

    case {'reconstruct'}
        path_data = load.path_to_data('_data');
        p = fullfile(path_data, 'Reconstruct', filesep);
        
    case 'stats'
        path_data = load.path_to_data('_data');
        p = fullfile(path_data, 'Stats', filesep);
        
        
    case {'impale_data', 'meas'}
        [before_root, ~, ~] = fileparts(rootpath);
        p = fullfile(before_root, '.data', filesep);
       
        
    case 'wav_spch_36sec'
        p_stimulus = load.path_to_data('stimulus');
        p = fullfile(p_stimulus, 'Spch_(36)sec');
        
    case 'timit_females_fs(100khz)'
        p_stimulus = load.path_to_data('wav_spch_36sec');
        p = fullfile(p_stimulus, 'TIMIT_Females_Fs(100kHz)');
        
    case 'timit_males_fs(100khz)'
        p_stimulus = load.path_to_data('wav_spch_36sec');
        p = fullfile(p_stimulus, 'TIMIT_Males_Fs(100kHz)');      
        
    case 'timit_females_fs(16khz)_shortlength'
        p_stimulus = load.path_to_data('wav_spch_36sec');
        p = fullfile(p_stimulus, 'TIMIT_Females_Fs(16kHz)_ShortLength');
        
    case 'timit_males_fs(16khz)_shortlength'
        p_stimulus = load.path_to_data('wav_spch_36sec');
        p = fullfile(p_stimulus, 'TIMIT_Males_Fs(16kHz)_ShortLength');           
        
        
    case 'raw'  % not on the HD, so it's "Hard-coded"
        %computer_name = getenv('computername');
        computer_name = computer;        
        if strcmpi('PCWIN64', computer_name)    %if strcmpi('OdedAlienware', computer_name)
            % OdedAlienware
            p = 'D:\Dropbox (Partners HealthCare)\DataBases\myMeas\Rabbits\!RAW';
        elseif strcmpi('MACI64', computer_name)  % elseif strcmpi('ODEDSDELL', computer_name)
            % Oded's external HD
            p = '/Users/ob993/Dropbox (Partners HealthCare)/DataBases/myMeas/Rabbits/!RAW';        
        else        
            % Apollo drive
            warning('--> Can''t find this COMPUTER! Trying to get the RAW path from APOLLO!');
            p = '//apollo/ent_archive/Delgutte_Archive/obz/Rabbits/C74/!RAW/';
        end                
        
    otherwise
        error('--> ERROR in [path_to_data.m]: Unrecognized DATA_TYPE option! (data_type = %s)',...
            data_type);
end
     

if ~isempty(verbose)
    fprintf('\n--> [path_to_data.m]:\n');
    fprintf('---> path  : %s\n', p);
    fprintf('---> is valid?: %d\n', ~isempty(dir(p)));
    fprintf('---> dir: \n');
    dir(p);
end

dir_p = dir(p);





