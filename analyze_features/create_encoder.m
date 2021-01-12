% create_encoder.m
%
%

clc
fignum = 11;
verbose = 1;

setup_environment('../');





%% Define the file information
finfo.type = 'MUA'   	% {'SU', 'MUA'}
% aux.cprintf('Keywords', '\n-> Data type: *** %s ***\n', finfo.type);
if strcmpi('MUA', finfo.type)    
    % *** DEFAULT ***
    % BINWIDTH : 5 ms
    % CAUSALITY: true
    % %{
    finfo.n_neurons  = 50; % [10 25 50 103 150 241];  	% # of units to load from existing files
    finfo.trainDRR   = '3';
    finfo.date       = '10-Jan-2021';        
    finfo.path       = '../_data/reconstruct/';
    fn_template      = 'reconstruct_%s_(%s)_units(%d)_bw(5)ms_algo(regression)_fbands(30)_splits(12)_lags(30)ms_cau(0)_trainDRR(%s).mat';
    %}
  
    
elseif strcmpi('SU', finfo.type)
    % BINWIDTH : 5 ms
    % CAUSALITY: true
    finfo.n_neurons  = [10 25 50 103];   	% # of units to load from existing files
    finfo.trainDRR   = '3';
    finfo.date       = '10-Jan-2021';        
    finfo.path       = '../_data/reconstruct/';
    fn_template       = 'reconstruct_%s_(%s)_units(%d)_bw(5)ms_algo(regression)_fbands(30)_splits(12)_lags(30)ms_cau(0)_trainDRR(%s).mat';

    
else
    error('-> Unrecognized measurement type (%s)!!',  finfo.type);
end
n_neuron_list = length(finfo.n_neurons);



% === Load first batch to get preliminary info
% finfo.fn = sprintf(fn_template, finfo.type, finfo.date, finfo.n_neurons(1), finfo.binwidth, ...
%     finfo.n_bands, finfo.win_size_ms, finfo.n_splits, finfo.trainDRR);
finfo.fn = sprintf(fn_template, finfo.type, finfo.date, finfo.n_neurons(1), finfo.trainDRR);
warning off
data     = load(fullfile(finfo.path, finfo.fn));
obj_list = data.obj_list;
spec_st  = data.spec_st;
warning on

% Make sure that the DRY index is the right one used in the loaded file!
%assert(drr.dry == obj_list{1,1}.train_drr);

% Sampling frequency along the time axis
binwidth    = spec_st.binwidth;     % (ms)
fs          = 1/(1e-3*binwidth);    % (Hz)
f           = spec_st.f;            % (Hz)
win_size_ms = spec_st.win_size_ms;  % (ms)

% OBJ_LIST
% * Dim1: # of DRR cases checked 
% * Dim2: # of  simulations 
% [n_drr, n_JK] = size(obj_list);
% assert(5 == n_drr);
% n_splits = n_JK;
n_bands  = obj_list{1}.n_bands;
% 




%% create the encoder
drr     = get_DRR_list_and_indices;
n_drr   = drr.n_drr;                  % # DRRs of used 
idx_dry = drr.ordered(1);
idx_drr = drr.ordered(5);

lags_ms     = 50;                       % (ms) the time duration of each patch 
binwidth    = spec_st.binwidth;         % (ms) 
n_lags      = fix(lags_ms/binwidth);    % (smp) number of samples of each patch (along the time duration)

Hdry = squeeze( data.H_units(:,idx_dry,:) );
Hdrr = squeeze( data.H_units(:,idx_drr,:) );

n_units = size(Hdry, 2);

n_splits = 100;
test_grp = 22;
[Htrn, Htst, ~, ~, split_st] = train_test_split(Hdry', randn(size(Hdry,1),1), ...
    'n_splits', n_splits, 'test_grp', test_grp);
Mtrn = im2col(Htrn, [n_units, n_lags], 'sliding');
% Mtst = im2col(Htst, [n_units, n_lags], 'sliding');
Mtrn = zca(Mtrn);

[~, Htst_] = train_test_split(Hdrr', randn(size(Hdrr,1),1), ...
    'n_splits', n_splits, 'test_grp', test_grp);
% Mtrn = im2col(Htrn, [n_units, n_lags], 'sliding');
Mtst_ = im2col(Htst_, [n_units, n_lags], 'sliding');
Mtst_ = zca(Mtst_);

Sdry = spec_st.Sft{idx_dry};
[Strn, Stst, ~, ~, split_st] = train_test_split(Sdry, randn(size(Sdry,2),1), ...
    'n_splits', n_splits, 'test_grp', test_grp);
Strn = im2col(Strn, [n_bands, n_lags], 'sliding');
Stst = im2col(Stst, [n_bands, n_lags], 'sliding');


Sdrr = spec_st.Sft{idx_drr};
[Strn_, Stst_, ~, ~, split_st] = train_test_split(Sdrr, randn(size(Sdrr,2),1), ...
    'n_splits', n_splits, 'test_grp', test_grp);
Strn_ = im2col(Strn_, [n_bands, n_lags], 'sliding');
Stst_ = im2col(Stst_, [n_bands, n_lags], 'sliding');



%%
k = 41
x = Mtst_(:, k);
enc = Mtrn' * x;
[~, idx_pp] = sort(enc, 'descend');     % posterior peak

figure(1);
plot(enc);
hold on
plot(idx_pp(1), enc(idx_pp(1)), 'ro');
hold off

Shat = zeros(n_bands*n_lags, 1);
for jj = 1:n_units
    Shat = Shat + Strn(:,idx_pp(jj));   
end


Rx = @(X) zca( reshape(X, n_bands, n_lags) );

figure(2);
clf;
subplot(1,3,1);
% imagesc([Rx(Stst_(:,k)), Rx(Shat), Rx(Stst(:,k))]);
imagesc( Rx(Stst_(:,k)) );
subplot(1,3,2);
imagesc( Rx(Shat) );
cc_ = corrcoef(Stst_(:,k), Shat);
cc = corrcoef(Stst(:,k), Shat);
title(sprintf('CC(est-drr): %.2f ----------  CC(est-dry): %.2f', cc_(1,2), cc(1,2)));
subplot(1,3,3);
imagesc( Rx(Stst(:,k)) );







