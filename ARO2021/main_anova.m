% main_anova.m
% 
% Description:
% Perform ANOVA analysis for the DRRs vs SUs & MUAs; which one is "better"?
%
% ref: https://www.mathworks.com/help/stats/two-way-anova.html
%


close all;
clc;
clear all;


fignum = 10;
verbose = 1;

setup_environment('../');

drr = get_DRR_list_and_indices;



%% Load RAW data
%   struct with fields: 
%       stats: [1×1 struct]
%     tbl_MUA: [241×21 table]
fn.path = load.path_to_data('Analysis');
fn.file = 'stats_raw(MUA)_(19-Feb-2021)_BW(5)ms_duration(36)sec_units(241).mat';
data.mua.raw = load( fullfile(fn.path, fn.file) );


% MUA
%   struct with fields:
% 
%         CCstats: [50×5×12 double]
%         CCunits: [50×50 double]
%        H_sorted: [7200×5×50 double]
%              fn: [1×1 struct]
%     sorted_list: [1×241 double]
%         spec_st: [1×1 struct]
%         stim_st: [1×1 struct]
%        tbl_data: [241×20 table]
fn.path = load.path_to_data('Analysis');
fn.file = 'CCstats_MUA_(14-Feb-2021)_units(50)_bw(5)ms_algo(regression)_fbands(30)_splits(12)_lags(30)ms_cau(0)_trainDRR(4)';
fn.file = [fn.file, '.mat'];
data.mua.drr4 = load( fullfile(fn.path, fn.file) );


fn.path = load.path_to_data('Analysis');
fn.file = 'CCstats_MUA_(13-Feb-2021)_units(50)_bw(5)ms_algo(regression)_fbands(30)_splits(12)_lags(30)ms_cau(0)_trainDRR(3)';
fn.file = [fn.file, '.mat'];
data.mua.dry = load( fullfile(fn.path, fn.file) );


% SU
fn.path = load.path_to_data('Analysis');
fn.file = 'CCstats_SU_(13-Feb-2021)_units(50)_bw(5)ms_algo(regression)_fbands(30)_splits(12)_lags(30)ms_cau(0)_trainDRR(3)';
fn.file = [fn.file, '.mat'];
data.su.dry = load( fullfile(fn.path, fn.file) );




%% dry \ drr4
figure(1);
subplot(1,2,1);
plot( std(squeeze(data.mua.dry.CCstats(:,:,1)),[],2), '.-' )

subplot(1,2,2);
ix = 49
plot(data.mua.dry.CCstats(ix,drr.ordered,1)')
hold on


%% drr4 \ drr4
figure(2); 
subplot(1,2,1);
plot( std(squeeze(data.mua.drr4.CCstats(:,:,1)),[],2), '.-' )

subplot(1,2,2);
ix = 1
plot(data.mua.drr4.CCstats(ix,drr.ordered,1)')
hold on






%%


%{
% Lets create a dummy data base:
n_rep = 30;     % number of repetitions
n_drr = drr.n_drr;
A = randn(n_rep, n_drr);   

bias = 1;
for k = 1:n_drr
    A(:,k) = A(:,k) + k*bias;    
end

figure(1); clf;
plot_dotbox(A);

[p, tbl, stats] = anova2(A, n_rep, 0);
tbl = cell2table(tbl);
tbl

figure(2); clf;
c = multcompare(stats);
c
%}



%% Compare between DRRs of MUAs
drr_1 = 3;      % DRR condirion #1
drr_2 = drr_1;      % DRR condirion #2 --- DRR = -8.2 dB
M1 = squeeze(data.su.dry.CCstats(:,drr_1,:))';     
M2 = squeeze(data.mua.dry.CCstats(:,drr_2,:))';     
n_rep = size(M1, 1);    % # of repetitions in each case

% Prepare Data for Balanced Two-Way ANOVA
mtx_anova2 = [M1(:), M2(:)];

[p, tbl, stats] = anova2(mtx_anova2, n_rep, 0);
tbl



%% Compare between DRRs of MUAs
drr_1 = 3;      % DRR condirion #1
drr_2 = 4;      % DRR condirion #2 --- DRR = -8.2 dB
M1 = squeeze(data.mua.drr4.CCstats(:,drr_1,:))';     
M2 = squeeze(data.mua.drr4.CCstats(:,drr_2,:))';     
n_rep = size(M1, 1);    % # of repetitions in each case

% Prepare Data for Balanced Two-Way ANOVA
mtx_anova2 = [M1(:), M2(:)];

[p, tbl, stats] = anova2(mtx_anova2, n_rep, 0);
tbl



%% Compare between DRRs of MUAs
drr_1 = 3;      % DRR condirion #1
drr_2 = 4;      % DRR condirion #2
M1 = squeeze(data.mua.dry.CCstats(:,drr_1,:))';
M2 = squeeze(data.mua.dry.CCstats(:,drr_2,:))';
n_rep = size(M1, 1);

% Prepare Data for Balanced Two-Way ANOVA
mtx_anova2 = [M1(:), M2(:)];

[p, tbl, stats] = anova2(mtx_anova2, n_rep, 0);
tbl



























