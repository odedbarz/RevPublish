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
    case 'reverberation'
        p = [rootpath, filesep];
        
    case 'meas'
        one_level_up_path = fileparts(rootpath);
        p = [one_level_up_path, filesep, '.data', filesep];
        
    case 'stimulus'
        p = [rootpath, filesep, '_data', filesep, 'Stimulus', filesep];             
        
    case 'wav_spch_36sec'
        p_stimulus = load.path_to_data('stimulus');
        p = [p_stimulus, 'Spch_(36)sec', filesep];
        
    case 'data'
        p = [rootpath, filesep, '_data', filesep];
        
    case 'impale_data'
        [codeOnCloud_path, ~, ~] = fileparts(rootpath);
        p = [codeOnCloud_path, filesep, '.data', filesep];
        
    case 'raw'  % not on the HD, so it's "Hard-coded"
        computer_name = getenv('computername');
        if strcmpi('OdedAlienware', computer_name)
            % OdedAlienware
            p = 'D:/DataBases/myMeas/Rabbits/!RAW/';

        elseif strcmpi('ODEDSDELL', computer_name)
            % Oded's external HD
            p = 'D:/oded_backups/myDataBases/DataBases 03-17-2020/myMeas/Rabbits/!RAW/';
        
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





