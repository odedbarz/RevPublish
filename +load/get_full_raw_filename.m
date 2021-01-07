function [raw_session, path_raw_session, path_root_raw] = ...
    get_full_raw_filename( fn_session, path_root_raw, verbose )

if 2 > nargin
    % path_root_raw = 'D:\DataBases\myMeas\Rabbits\!RAW\';  % odedalienware
    path_root_raw = 'D:\oded_backups\myDataBases\DataBases 03-17-2020\myMeas\Rabbits\!RAW\';    % odedlaptop   
end

if 3 > nargin
    verbose = [];
end


% Get just the filename, without the path & extention, if provided
[~, fn_session, ~] = fileparts(fn_session);

% Cut out the session #, Unit #, and measurement #
str_idx = regexp(fn_session, '-');
if ~isempty(str_idx)
    fn_session_name = fn_session(1:(str_idx(1)-1));
end

% Path to the raw session 
path_raw_session = [path_root_raw, 'et_', fn_session_name, filesep];

% Filename of the desired raw session
raw_session = [path_raw_session, fn_session, '.0.et'];

if ~isempty(verbose)
    fprintf('\n--> [path_raw_session.m]:\n');
    fprintf('---> path  : %s\n', path_raw_session);
    fprintf('---> isdir?: %d\n', isfolder(path_raw_session));
    fprintf('---> dir: \n');
    dir(path_raw_session);
    fprintf('--->\n');
    dir(raw_session);
end


