% main_all_meas_stats.m
%

close all;
clc;
clear all;

fignum = 10;
verbose = 1;

setup_environment('../');
debug_flag = 0;


%% Load data
data_type   = 'MUA';       % {'SU', MUA'}
fn.load.path= load.path_to_data('_data');
data_type   = upper(data_type);

fn.load.file_template = 'data_%s_(08-Nov-2021)_bw(5)_fbands(30)_win(NaN)ms_spec(gammatone)_ALL-MEAS';
fn.load.file = sprintf(fn.load.file_template, data_type);
fn.load.fullfile = fullfile( fn.load.path, [fn.load.file, '.mat'] );
dummy = load(fn.load.fullfile, 'tbl_MUA_all');
tbl_MUA_all = dummy.tbl_MUA_all;



%% Get a list of SU units
fn.load.file = sprintf(fn.load.file_template, 'SU');
fn.load.fullfile = fullfile( fn.load.path, fn.load.file );
data   = load(fn.load.fullfile, 'tbl_SU');
neurons= data.tbl_SU.neuron;



%%
tbl_MUA_all = tbl_MUA_all(neurons, :);
