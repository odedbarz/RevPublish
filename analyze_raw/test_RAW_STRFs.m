%
% test_RAW_STRFs.m
%
%
% Description:
% Learns the STRFs of the raw data from all recorded channels.
% 

% Add the path of the StimViewerGUI GUI
%addpath(genpath('../../StimViewerGUI'));
addpath('../');

FigSetup;

fignum = 11;


%% Get the stimuli & spectrograms
if ~exist('tbl_all', 'var') || isempty(tbl_all)
    data_st.path	= load.path_to_data; 
    data_st.fn      = 'LabNoteBook_Reverb.xlsx';
    data_st.operator= 'OBZ';
    data_st.ID      = '';
    data_st.measType= 'Spch';   % {'Spch', 'ns_Spch_fc1kHz', 'ns_Spch_fc4kHz'}
    data_st.session = [];     
    data_st.unit    = [];      
    data_st.measNum = [];     
    
    tbl_all = load.data_table(data_st);
end

%{
suffix      = 'Spch';   % S.measParam.Suffix;
duration_sec= 36;       % S.stimChans{1}.Source.numTokens;
fn          = get_stimulus_info(suffix, duration_sec);
Fs          = 100e3;    % (Hz)
stim = load.stimuli(Fs, fn.path, fn.template);

% (cell --> matrix)
stim.Y = [stim.Y{:}];
%}

% Make the BINWIDTH corresponds to the MUA & LFP (sampling rate of 500 Hz)
binwidth = 2;   % (ms)

%
[stim, spec, ~, data_st]  = load.stimulus_and_spectrogram(tbl_all, 1,...
    'binwidth', binwidth);
fprintf('--> binwidth: %g ms\n', spec.binwidth);

% % Loads PSTHs
% S_list = load.response( tbl_all(8,:) );
% psth = calc_PSTHs(S_list, spec.binwidth);


%% Impale files (MAT files)

% OdedAlienware
% path_root_mat = 'D:\DataBases\myMeas\Rabbits\OBZ_C74_022\';

% Oded's external HD
path_root_mat = 'C:\Users\barzelo\Google Drive\codeOnCloud\.data\';


%% Raw files (0.et files)
% OdedAlienware
% path_root_raw = 'D:\DataBases\myMeas\Rabbits\!RAW\';

% Oded's external HD
path_root_raw = 'D:\oded_backups\myDataBases\DataBases 03-17-2020\myMeas\Rabbits\!RAW';

% Apollo drive
% path_root_raw = '\\apollo\ent_archive\Delgutte_Archive\obz\Rabbits\C74\!RAW\et_OBZ_C74_020\';


%% Load Impale's structure
clear fn
fn.session     =  'OBZ_C74_019';                        % 'OBZ_E16_002'; 'OBZ_C74_029';
fn.file        = [fn.session, '-3-7-Spch_MERGED'];   	% '-2-2-Spch'  ; '-2-2-Spch'
seps           = regexp(fn.file, '-');
fn.session_dir = fn.file(1:seps(1)-1);

% Impale
fn.path.Impale = [path_root_mat, fn.session_dir, filesep];
fn.Impale      = [fn.path.Impale, fn.file, '.mat'];
assert(~isempty(dir(fn.Impale)), sprintf('--> Can''t find %s', fn.Impale));

% Impale structure
warning off
S = load(fn.Impale);
warning on
% S = load.Impale_struct( fn.Impale );


% Stimulus duration
assert( 1e-3*stim.duration_ms == S.stimChans{1}.Source.numTokens );
duration_ms  = stim.duration_ms;
duration_sec = 1e-3*duration_ms;



%% Load the RAW file
fn_raw_tbl = [path_root_raw, filesep, S.info.rawDataDir, filesep, S.info.rawFileName];
assert(~isempty(dir(fn_raw_tbl)), '--> ERROR: can''t find this file!!');

raw_st = load(fn_raw_tbl, '-mat');
assert(duration_ms == raw_st.rawdata.t_duration);



%%
fignum = 11;

inner     = 1;
outer     = 1;
channel   = 5;    %1+fix(spikechan/6);     % channels number
spikechan = 1 + (channel-1)*6;
%trial     = [];      % trial number to look at

tbl_raw = raw_st.tbl;

idx_tbl = (tbl_raw.outer == outer) &...
    (tbl_raw.inner == inner) &...
    (tbl_raw.channel==channel);
assert(1<=nnz(idx_tbl), '--> Did not find such measurement!');

n_trials = nnz(idx_tbl);

len_x = length(tbl_raw.x{1});    % (samples)


% Time axis
%t = linspace(0, duration_sec, length(tbl_raw.x{1}))';   % (sec) 

% Table to one matrix
X = [tbl_raw.x{idx_tbl}];

%  Multiunit activity
[mua, pars_mua] = MUA(X, raw_st.sr); 
mua_avg = mean(mua,2);

% The response to use for the STRF 
response = mua_avg;

% Time axis
t = linspace(0, duration_sec, length(response))';   % (sec) 

%

% Split the data into train & test sets
% A division of the concatenated stimulus into the TIMIT WAVs
len_resp = length(response);
n_groups = 12;
assert( fix(len_resp/n_groups) == len_resp/n_groups, '--> Please choose number of groups that divide the number of samples!' );
n_sample_in_grp = fix(len_resp/n_groups);
groupIndex = cumsum((0 == mod((1:length(response))-1, n_sample_in_grp))) ;

%     '@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@'
%     a = find(groupIndex ==12);
%     b = find(groupIndex ~=12);
%     groupIndex = groupIndex([b, a(1:1200)]);

% Debug: 
%plot(groupIndex);
% % make sure that all group are of equal size
%arrayfun(@(ii) nnz(groupIndex==ii), 1:n_groups)     


% Revb.          20%        80%   100%        20%        80%     100%           
% Dist.         1.5m       1.5m   1.5m       3.0m       3.0m     3.0m
drr_list = {'-2.5 dB',  '9.4 dB', 'Dry', '-8.2 dB',  '4.8 dB'};
dist_list =[     1.5,       1.5,   1.5,       3.0,       3.0];
revb_list =[      20,        80,   100,        20,        80];

idx.use  = inner * outer;  
aux.cprintf('g', '--> Sft used: *** %s ***\n', drr_list{idx.use});

%
Sft = spec.Sft{idx.use};
%psth= psth(:,idx.use);

% Split the data into train & test sets
clear data
train.grp  = 1:(n_groups-1);
train.idx  = strfpkg.get_selected_chunk_indices(train.grp, groupIndex);
train.Sft  = Sft(:, train.idx);
%train.resp = psth(train.idx)';
train.resp = response(train.idx)';

test.grp  = n_groups;
test.idx  = strfpkg.get_selected_chunk_indices(test.grp, groupIndex);
test.Sft  = Sft(:, test.idx);
%test.resp = psth(test.idx);
test.resp = response(test.idx);

%{
clear data
data.groupIndex = groupIndex;
data.train.grp  = 1:11;

data.test.grp  = 12;
data.test.idx  = strfpkg.get_selected_chunk_indices(data.test.grp, data.groupIndex);
data.test.Sft  = wholeStim(data.test.idx,:)';
data.test.resp = wholeResponse(data.test.idx)';

% Print some information about the data
fprintf('\n');
fprintf('--> Total stimulus length: %g sec\n', sum(stimInfo.stimLengths));
%}



%% STRF
% (samples) the window to take for the xcorr (maxlags)
n_win = 2*39+1;     

% Ridge regression constant(s)
tol = 0.0051;

% calculate the strf(s)
[strf, strf_st] = strfpkg.strf(Sft, response, n_win,...
    'tol', tol,...
    'faxis', spec.f, ...
    'binwidth', spec.binwidth, ...     
    'fignum', []);


figure(fignum);     
fignum = fignum + 2;
clf;
strfpkg.plot_strf(strf, spec.f, binwidth);


% Estimate the response
r_est = strfpkg.reconstruct_response( test.Sft, strf,...
    'avg', strf_st.avg.psth, ...
    'normalize', max(test.resp));


figure(fignum);     
fignum = fignum + 2;
clf;
subplot(2,1,1);
plot(t(1:length(test.resp)), [test.resp, r_est]);

subplot(2,1,2);
zsc = @(x) (x-mean(x))./std(x);
plot(t(1:length(test.resp)), zsc([test.resp, r_est]));
ylabel('Z-score');
CC = corrcoef(test.resp, r_est);
title(sprintf('CC: %.3f', CC(1,2)));


