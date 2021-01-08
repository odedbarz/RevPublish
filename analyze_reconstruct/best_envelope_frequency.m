% best_envelope_frequency.m
% 
% Description:
% This script finds the best frequency (BF) according to CC between measurements 
% response and the envelope of the stimulus.


clc
fignum = 11;
verbose = 1;

setup_environment('../');


drr   = get_DRR_list_and_indices;
n_drr = drr.n_drr;                  % # DRRs of used 




%% Load data
%   Run [main_aggregate_MUA_data.m] again to update this file if needed
% 
data_type    = 'MUA';       % {'SU', MUA'}
data_type    = upper(data_type);
fn.load.path = '../_data';

switch data_type
    case 'SU'
        % Loads a struct with fields:
        %               H: [3600×6×150 double]
        %          S_list: {1×150 cell}
        %     neuron_list: [150×1 double]
        %         spec_st: [1×1 struct]
        %      tbl_impale: [437×20 table]
        fn.load.file = 'data_SU_(07-Jan-2021)_bw(5)_fbands(30)_win(NaN)ms_spec(gammatone).mat';       
        
    case 'MUA'
        %Loads a struct with fields:
        %               H: [7200×6×356 double]
        %        H_labels: [356×6 double]
        %     neuron_list: [356×1 double]
        %         spec_st: [1×1 struct]
        %         stim_st: [1×1 struct]
        %      tbl_impale: [437×20 table]        
        fn.load.file = 'data_MUA_(07-Jan-2021)_bw(5)_fbands(30)_win(NaN)ms_spec(gammatone).mat';        
        
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



% %% Get the valid measurements\columns
% % Desired stimulus duration to use
% duration_sec = 36;      % (sec) 
% 
% % Make sure that the loaded data is of the right duration
% assert(duration_sec == spec_st.duration_ms * 1e-3, '--> ERROR: You are using the wrong stimulus duration!');
% 
% % Select all measurement with the desired duration
% neuron_duration_list = tbl_impale.duration_sec == duration_sec;
% 
% % Select measurements that has all recorded sessions
% slc.valid_neuron_idx = n_drr <= sum( ~isnan( squeeze(sum(data.H,1)) ), 1)';
% slc.valid_neuron_idx = slc.valid_neuron_idx(:);
% 
% % Indices of both boolean conditions
% slc.valid_neuron_idx = slc.valid_neuron_idx & neuron_duration_list(data.neuron_list);
% 
% switch data_type
%     case 'SU'
%         % Make sure that all SUs are valid!
%         assert(all(1 == tbl_impale.SPK( data.neuron_list(slc.valid_neuron_idx) )),...
%             '--> ERROR: some of these n_units don''t contain a SU!');        
%         n_units = nnz(slc.valid_neuron_idx);
%         
%     case 'MUA'
%         n_units = nnz(slc.valid_neuron_idx);
%         
%     otherwise
%         error('--> Unrecognized DATA_TYPE!');
%         
% end
%         
% % A list of all "valid" n_units to use
% slc.valid_neurons = data.neuron_list(slc.valid_neuron_idx);
% 
% 
% % Get a valid list of neurons (duration & all-sessions available)
% H_valid = data.H(:,:,slc.valid_neuron_idx);
% assert( size(H_valid,3) >= max(n_units), ...
%     '--> You are asking for more neurons than are available in the dataset!');



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
    ydry =  squeeze(data.H(:,dry_idx,n));
    
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
tbl_BF = [table(tbl_data.Row, 'VariableNames', {'Row'}), table(BF), table(BF_cc), table(BF_pv)];

save([fn.load.fullfile(1:end-4), '_BF'], 'tbl_BF', 'spec_st', 'stim_st');
%}










