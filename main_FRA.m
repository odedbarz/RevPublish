%
% main_FRA.m.m
% 
% Description:
%   Reads all FRAs from the raw data.
% 
%

% viewer.plot_FRA(S);


clc
fignum = 11;
verbose = 1;

% setup_environment('../');
data_path = load.path_to_data('raw');


%% Load data
%   Run [main_aggregate_MUA_data.m] again to update this file if needed
% 
data_type   = 'SU';       % {'SU', MUA'}
fn.load.path= load.path_to_data('_data');
data_type   = upper(data_type);

fn.load.file_template = 'data_%s_(08-Jan-2021)_bw(5)_fbands(30)_win(NaN)ms_spec(gammatone).mat';
% fn.load.file_template = 'data_%s_(01-Nov-2021)_bw(1)_fbands(30)_win(NaN)ms_spec(gammatone).mat';
% fn.load.file_template = 'data_%s_(08-Nov-2021)_bw(5)_fbands(30)_win(NaN)ms_spec(gammatone).mat';
% fn.load.file_template = 'data_%s_(10-Dec-2021)_bw(100)_fbands(30)_win(NaN)ms_spec(gammatone).mat'

fn.load.file = sprintf(fn.load.file_template, data_type);
fn.load.fullfile = fullfile( fn.load.path, fn.load.file );
dummy = load(fn.load.fullfile, 'tbl_SU');
tbl_SU = dummy.tbl_SU;



%%
