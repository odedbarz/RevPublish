function [X_, rawdata] = raw_spch_wave( S, fn_raw ) 
%
%   function [X_scaled, rawdata] = raw_spch_wave( S, fn_raw )  
%
% Input:
%   mfull_filename: The raw filename (full path) to load.
%   metadata: structure with flags and fields about how manipulate the
%             loaded data.
%
% Output:
%   o: (struct) structure of the raw recorded data
%       y: (sum of (inner x outer x rep)) raw vectors
%       idx: the location (index) in the (inner x outer x rep) matrix
%       Fs: (Hz) Sampling rate of the output raw data
%
%
% Impale's data:
%   params: raw info parameters:
%       params.Fs         : Sampling rate (Hz)
%       params.scaleFactor: 
%       params.npts_bad   :
%       params.n_spikechan      : # of channels (nchan)
%
%
% Description:
%   The original file is taken from Brian's Raw_waveform_viewer GUI. 
%   
%
%

% Todo: get it out as a parameter in the class
et_header_bytes = 24;       % Set the header of the 0.et file

% How much bytes is an int16?
dummy = int16(1);
dummy_info = whos('dummy');
int16bytes = dummy_info.bytes;


%% Open the raw data file (.0.et)
assert(~isempty(dir(fn_raw)), sprintf('--> Can''t find this file:\n\t----> %s!\n', fn_raw));
fid = fopen(fn_raw, 'rb');
if fid < 0, error('-> [+load.raw_spch_wave.m]: Cannot open file for reading!'); end


%% Read meta-data from the et file's header
% * There are 4 numbers BEFORE the actual data starts
Fs          = fread(fid,1,'double');        % (Hz) Sampling rate;
scaleFactor = fread(fid,1,'double');
npts_read   = fread(fid,1,'int32');         % (smp) num of samples per trial
n_spikechan = fread(fid,1,'int32');         % (smp) num of recording channels; can be 3 (for 1 channel, the first two are identical and the third is zero) or 4 (for 4 channels);

% (Hz) Correct the sampling rate if needed (Usually I'm using 100k Hz)
% Fs_out = S.dacInt.SampleRate;  % (Hz) output sampling rate
if ischar(Fs), Fs = str2double(Fs); end
assert(Fs == S.dacInt.SampleRate, '--> ERROR: mismatched sampling frequencies!!!');


%% Time durations:
t_token = S.dacInt.TotalDuration;      % (ms) duration of each token
if ischar(t_token), t_token = str2double(t_token); end

% # of repeats for each token; e.g., 36 tokens (=n_tokens), each one is 
% 1000 ms (=t_token) in length, for a total speech stimulus duration of 36 sec 
n_tokens = S.stimChans{2}.Source.numTokens;     

% # of samples in each token; e.g., each token is 1000 ms long and the
% sampling rate in Impale is 100kHz (=Fs), so there are 100,000 samples in
% each token
smp_token = fix( units.ms2sec(t_token)*Fs );    % (smp) 


% Duration of one speech stimulus
t_duration = t_token .* n_tokens;   % (ms) 

% # of samples for each speech interval
smp_duration = fix( units.ms2sec(t_duration)*Fs );    % (smp)             

% * File length is the total length minus the 24 bytes from above. 
% * The number of trials is that number divided by nchan, npts, and 2 (#bytes in int16)
fseek(fid,0,'eof');         % get to end of file
filesize = ftell(fid);      % [bytes] get the total file size
datasize = (filesize - et_header_bytes)/int16bytes;  % [bytes] get the size of the data (without the header)

% # of trials per speech interval
n_trials = datasize/(n_spikechan*smp_duration);   % number of trials
assert(n_trials == fix(n_trials), '--> number of trials MUST be an integer!')

% # of repetirions for each token
n_repits_per_token = S.measParam.stopVal;

% # of trials per x1 full length speech stimulus  
n_trials_per_stimulus = n_repits_per_token/n_tokens; 

assert(0 == mod(n_trials,n_trials_per_stimulus),...
    '--> These two MUST share a common factor!');



%% Load the data with the PST data
fseek( fid, et_header_bytes, 'bof' );   % move to the beginning of the file + et_header_bytes
[X, n_loaded] = fread(fid, datasize, 'int16');
assert(n_loaded == datasize, '--> FREAD didn''t read all requested data!');
fclose(fid);

assert(datasize == n_loaded, 'Number of loaded samples isn''t equal to the number of requested samples!!');   

% Recale the measurements; X & X_scale are one lone vector that contains
% all tokens\trials
X_ = X/scaleFactor;

% For ONE channel, clean the reduncdant data; There are 3 same 
% repetitions; just one is needed
if 3 == n_spikechan     
    % Remove two repetitions out of the 3
    %X_ = load.clean_one_channel_redundancy( X_scaled, n_spikechan, smp_duration );
    X_ = reshape(X_, smp_token, []);
    X_ = X_(:,1:n_spikechan:end);
    assert(size(X_,2) == n_trials*n_tokens)
    
    % So, now X_ looks like this:
    %   X_: [1:smp_duration, 1:smp_duration, ..., 1:smp_duration]
    % with N_TRIALS repetitions
    X_ = X_(:);
        
    n_spikechan = 1;
end


%% More information about the recorded measurement
idx.inner = [S.SCL(:).innerIndex];
idx.valid = idx.inner > 0;
idx.inner = idx.inner(idx.valid);
n_inner   = max(idx.inner);
assert(n_inner >= 1);

idx.outer = [S.SCL(:).outerIndex];
idx.outer = max(1, idx.outer(idx.valid));
n_outer   = max(idx.outer);
assert(n_outer >= 1);

%idx.repIndex = [S.SCL(:).repIndex];



%% WAVE POLARITY: 
% set the polarity of all channels to have the same sign
polarity = [S.discrimSettings.settings(:).polarity];    % 0 or 1
polarity = (2*polarity - 1);    % length(polarity) == n_spikechan

%% # of channels
% Get only the ones with spikes
%syncchan = 0;
%spikechan_all = spiketools.get_all_spikechan(S, syncchan);

% Get all of them!
spikechan_all = channel_2_spikechan( 1:length(polarity) );


%% Set the output
rawdata.fn_raw = fn_raw;            % (str) raw filename   
rawdata.outerSeq = idx.outer; % S.outerSeq.master.values;
rawdata.outer_iorder = S.outerSeq.master.iorder;
rawdata.outerIndex = idx.outer;
rawdata.innerSeq = idx.inner; % S.innerSeq.master.values;
rawdata.inner_iorder = S.innerSeq.master.iorder;
rawdata.innerIndex = idx.inner;
rawdata.n_inner  = n_inner;
rawdata.n_outer  = n_outer;
rawdata.Fs = Fs;                    % (Hz) (1x1) sampling rate 
rawdata.t_token = t_token;          % (ms) duration of each token (S.dacInt.TotalDuration)
rawdata.n_tokens = n_tokens;        % # of repeats for each token  (S.stimChans{2}.Source.numTokens) 
rawdata.smp_token = smp_token;      % # of samples in each token
rawdata.t_duration = t_duration;    % (ms) Duration of one speech stimulus (t_token .* n_tokens)
rawdata.smp_duration = smp_duration;% # of samples of each speech interval
rawdata.n_trials = n_trials;        % # of trials
rawdata.n_trials_per_stimulus = n_trials_per_stimulus;
rawdata.spikechan_all = spikechan_all;
rawdata.n_spikechan = n_spikechan;	% # of recorded channels 
rawdata.polarity = polarity;
rawdata.S = S;                      % (struct) impale's structure
rawdata.datasize = datasize;        % (1x1) # of loaded data elements 


