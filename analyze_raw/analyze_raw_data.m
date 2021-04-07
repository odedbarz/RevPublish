%
% analyze_raw_data.m
%
% Description:
% Analyze the raw data file,
%   stats_raw(MUA)_(19-Feb-2021)_BW(5)ms_duration(36)sec_units(241).mat
%
% To create this file, analyze_raw/run main_raw_MUA_stats.m
%


clc
fignum = 10;
verbose = 1;

setup_environment('../');

drr   = get_DRR_list_and_indices;
n_drr = drr.n_drr;                  % # DRRs of used 
idx_dry = drr.ordered(1);


%% Load RAW data
%   struct with fields: 
%       stats: [1×1 struct]
%     tbl_MUA: [241×21 table]
fn.path = load.path_to_data('Analysis');
fn.file = 'stats_raw(MUA)_(19-Feb-2021)_BW(5)ms_duration(36)sec_units(241).mat';
load( fullfile(fn.path, fn.file) );



%%





