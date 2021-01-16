% create_BF_table.m
% 
% Description:
% This script finds the best frequency (BF) according to CC between measurements 
% response and the envelope of the stimulus.


clc
fignum = 11;
verbose = 1;

setup_environment('../');



%% Load data
%   Run [main_aggregate_MUA_data.m] again to update this file if needed
% 
data_type    = 'SU';       % {'SU', MUA'}
data_type    = upper(data_type);
fn.load.path = '../_data';

drr   = get_DRR_list_and_indices;
n_drr = drr.n_drr;                  % # DRRs of used 

switch data_type
    case 'SU'
        % Loads a struct with fields:
        %               H: [3600×6×150 double]
        %          S_list: {1×150 cell}
        %     neuron_list: [150×1 double]
        %         spec_st: [1×1 struct]
        %      tbl_impale: [437×20 table]
        fn.load.file = 'data_SU_(08-Jan-2021)_bw(5)_fbands(30)_win(NaN)ms_spec(gammatone).mat';       
        
    case 'MUA'
        %Loads a struct with fields:
        %               H: [7200×6×356 double]
        %        H_labels: [356×6 double]
        %     neuron_list: [356×1 double]
        %         spec_st: [1×1 struct]
        %         stim_st: [1×1 struct]
        %      tbl_impale: [437×20 table]        
        fn.load.file = 'data_MUA_(08-Jan-2021)_bw(5)_fbands(30)_win(NaN)ms_spec(gammatone).mat';        
        
    otherwise
        error('--> Unrecognized DATA_TYPE!');
        
end
fn.load.fullfile = fullfile( fn.load.path, fn.load.file );
data    = load(fn.load.fullfile);
spec_st = data.spec_st;
stim_st = data.stim_st;
tbl_data= data.(sprintf('tbl_%s', data_type));
n_units = height(tbl_data);
duration_sec = 36;      % (sec) 


aux.vprint(verbose, '--> [main_loopover_units.m] Loading file:\n\t...<%s>\n', fn.load.file);
aux.vprint(verbose, '-> data_type: %s\n', data_type);



%%
binwidth       = spec_st.binwidth;     % (ms)
n_bands        = spec_st.n_bands;
win_size_ms    = spec_st.win_size_ms; 


if verbose
    aux.cprintf('UnterminatedStrings', '\n    Data:\n');
    aux.cprintf('UnterminatedStrings', '--> data_type   : %s\n', data_type);
    aux.cprintf('UnterminatedStrings', '--> duration_sec: %g ms\n', duration_sec);
    aux.cprintf('UnterminatedStrings', '    Reconstruction:\n');
    aux.cprintf('UnterminatedStrings', '--> binwidth    : %g ms\n', binwidth);
    aux.cprintf('UnterminatedStrings', '--> n_bands     : %g\n', n_bands);
    aux.cprintf('UnterminatedStrings', '--> win_size_ms : %g ms\n', win_size_ms);
end

    
    
%% Loop over all n_units
dry_idx = drr.ordered(1);
Sdry    = spec_st.Sft{dry_idx};     % fry spectrogram
n_freq  = size(Sdry, 1); 

% Best-frequency correlation coefficient
BF    = nan(n_units, 1);
BF_cc = nan(n_units, 1);
BF_pv = nan(n_units, 1);
alpha = 0.05;   % significance

for n = 1:n_units    
    % extract only the dry measurements 
    ydry =  squeeze( data.H(:,dry_idx,n) );
    
    Rbest = -1;         % best correlation
    Pbest = nan;        % best p-value    
    kbest = nan;
    for k = 1:n_freq
        [Rn, Pn] = corrcoef(Sdry(k,:), ydry, 'alpha', alpha );
        
        if Rbest < Rn(1,2)
            Rbest = Rn(1,2);
            Pbest = Pn(1,2);
            kbest = k;
        elseif Rbest == Rn(1,2) && Pbest > Pn(1,2)
            Pbest = Pn(1,2);
            kbest = k;
        end
        
        
    end
    
    % Get the best frequency over all envelopes of the spectrogram
    BF(n)    = spec_st.f(kbest);    % (Hz)
    BF_cc(n) = Rn(1,2);             % correlation coefficient of the best frequency
    BF_pv(n) = Rn(1,2);             % p-value of the correlation coefficient
    
end


%% Add BF as columns into the measurement table
% Create a new table with all the information of the measurements and save it
neuron = tbl_data.neuron;
tbl_BF = [table(neuron), table(BF), table(BF_cc), table(BF_pv)];

save([fn.load.fullfile(1:end-4), '_BF'], 'tbl_BF', 'spec_st', 'stim_st');
%}










