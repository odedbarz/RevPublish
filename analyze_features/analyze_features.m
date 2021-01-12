% analyze_features.m
% 
%

clc
fignum = 11;
verbose = 1;

setup_environment('../');



%% Load spectrogram data
%{
data_type = 'MUA';       % {'SU', MUA'}
path_stim = '../_data';
data_type = upper(data_type);
switch data_type
    case 'SU'
        % Loads a struct with fields:
        %               H: [3600×6×150 double]
        %          S_list: {1×150 cell}
        %     neuron_list: [150×1 double]
        %         spec_st: [1×1 struct]
        %      tbl_impale: [437×20 table]
        fn_data = 'data_SU_(08-Jan-2021)_bw(5)_fbands(30)_win(NaN)ms_spec(gammatone).mat';       
        
    case 'MUA'
        %Loads a struct with fields:
        %               H: [7200×6×356 double]
        %        H_labels: [356×6 double]
        %     neuron_list: [356×1 double]
        %         spec_st: [1×1 struct]
        %         stim_st: [1×1 struct]
        %      tbl_impale: [437×20 table]        
        fn_data = 'data_MUA_(08-Jan-2021)_bw(5)_fbands(30)_win(NaN)ms_spec(gammatone).mat';        
        
    otherwise
        error('--> Unrecognized DATA_TYPE!');
        
end
fn_2_load = fullfile( path_stim, fn_data );
dummy = load(fn_2_load);
spec_st = dummy.spec_st;
n_bands = spec_st.n_bands;
%}



 %% Define the file information
finfo.type = 'SU'   	% {'SU', 'MUA'}
aux.cprintf('Keywords', '\n-> Data type: *** %s ***\n', finfo.type);

% trainDRR = [3, 5];
% train_DRR_labels = [];
% train_DRR_idx = [];
% for k = 1:length(trainDRR)
%     train_DRR_labels = [train_DRR_labels, drr.labels{trainDRR(k)}, ', '];
%     train_DRR_idx    = [train_DRR_idx, num2str(trainDRR(k)), ', '];
% end
% train_DRR_labels = train_DRR_labels(1:end-2);
% train_DRR_idx = train_DRR_idx(1:end-2);
% aux.cprintf('Keywords', '-> trainDRR: [%s] (%s)\n\n', train_DRR_labels, train_DRR_idx);

if strcmpi('MUA', finfo.type)    
    % *** DEFAULT ***
    % BINWIDTH : 5 ms
    % CAUSALITY: true
    % %{
    finfo.n_neurons  = [10 25 50 103 150 241]; %[1 5 10 25 50 100 150] % 241];  	% # of units to load from existing files
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
dummy    = load(fullfile(finfo.path, finfo.fn), 'obj_list', 'spec_st');
obj_list = dummy.obj_list;
spec_st  = dummy.spec_st;
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
[n_drr, n_JK] = size(obj_list);
assert(5 == n_drr);
n_splits = n_JK;
n_bands  = obj_list{1}.n_bands;


if verbose
    aux.cprintf('UnterminatedStrings', '--> binwidth   : %g ms\n', binwidth);
    aux.cprintf('UnterminatedStrings', '--> n_bands    : %g\n', n_bands);
    aux.cprintf('UnterminatedStrings', '--> win_size_ms: %g ms\n', win_size_ms);
    aux.cprintf('UnterminatedStrings', '--> n_splits   : %d\n', n_splits);
end




%% Load the dictioary
path_2_dict = fullfile( load.path_to_data('data'), 'Dictionary');
fn_dic  = 'DIC_MUA_(11-Jan-2021)_type(NNMF)_lags(30)ms_bw(5)_fbands(30)_win(103)ms_spec(ammatone).mat';
fn_2_load = fullfile( path_2_dict, fn_dic);
dummy   = load(fn_2_load);    % --> 'D', 'dict_st', 'RS', 'spec_st'
D       = dummy.D;
RS      = dummy.RS;
n_lags  = dummy.dict_st.input_pars.n_lags;
dict_st = dummy.dict_st;
n_samples = size(RS, 2);
n_features = size(RS, 1);

fprintf('\nLoading dictionary:\n');
fprintf('-> filename: %s\n', fn_dic);
fprintf('-> path    : %s\n', path_2_dict);



%%
drr     = get_DRR_list_and_indices;
n_drr   = drr.n_drr;                  % # DRRs of used 
idx_dry = drr.ordered(1);

idx_test = 4;

% Sdry = spec_st.Sft{idx_dry};
% Xdry = im2col(Sdry, [n_bands, n_lags], 'sliding');

% Stest = spec_st.Sft{idx_test};
% Stest = Stest(:, 1:700);
% Xtest = im2col(Stest, [n_bands, n_lags], 'sliding');

n_speaker = 2;
Xdry  = obj_list{n_speaker, drr.ordered(1)}.X_est;      '### DEBUG ###'
Xtest = obj_list{n_speaker, drr.ordered(5)}.X_est;      '### DEBUG ###'

n_tests = 50;
encoding = nan(n_features, n_tests);
enc_score = nan(size(RS,2), n_tests);


rng('default'),     '### DEBUG ###'
s = rng       ,     '### DEBUG ###'

CCdry = zeros(1, n_tests);
CC = zeros(1, n_tests);
idx = randperm(n_tests);
for k = 1:n_tests
    pk = Xtest(:,idx(k)+[1:n_lags]-1);   % the k'th patch 
    pk = pk(:);
    
    encoding(:,k) = D'*pk/length(pk);
    %enc_score(:,k) = sum(abs(encoding(:,k) - RS).^2)';
    
    dummy = corrcoef(pk, D*encoding(:,k));
    CC(k) = dummy(1,2);
    
    % dry
    pk_dry = Xdry(:, idx(k)+[1:n_lags]-1);
    dummy = corrcoef(pk_dry, D*encoding(:,k));
    CCdry(k) = dummy(1,2);
end
imagesc(enc_score);
xlabel('Num. of Tests');
ylabel('Features');



%%
'### DEBUG ###'
a = reshape(pk, 30, n_lags);
adry = reshape(pk_dry, 30, n_lags);
b = reshape(D*encoding(:,end), 30, n_lags);
imagesc( zca([a, b, adry]) );
cc = corrcoef(a(:), b(:));
title( sprintf('CC: %.2f', cc(1,2) ) );






















