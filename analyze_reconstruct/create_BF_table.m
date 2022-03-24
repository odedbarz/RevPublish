% create_BF_table.m
% 
% Description:
% This script finds the best frequency (BF) according to CC between measurements 
% response and the envelope of the stimulus.


clc
fignum = 11;
verbose = 1;

setup_environment('../');

save_results = true


%% Load data
%   Run [main_aggregate_MUA_data.m] again to update this file if needed
% 
data_type   = 'MUA';       % {'SU', MUA'}
fn.load.path= load.path_to_data('_data');
data_type   = upper(data_type);

fn.load.file_template = 'data_%s_(13-Jan-2022)_bw(5)_fbands(30)_win(NaN)ms_spec(gammatone-only)';
% fn.load.file_template = 'data_%s_(13-Jan-2022)_bw(5)_fbands(30)_win(NaN)ms_spec(gammatone-SyncFilter)';
% fn.load.file_template = 'data_%s_(13-Jan-2022)_bw(5)_fbands(30)_win(NaN)ms_spec(meddis)';
% fn.load.file_template = 'data_%s_(13-Jan-2022)_bw(5)_fbands(30)_win(NaN)ms_spec(carney)';
% fn.load.file_template = 'data_%s_(25-Jan-2022)_bw(1)_fbands(90)_win(NaN)ms_spec(gammatone)';
%                        fn.load.file_template = 'data_%s_(08-Jan-2021)_bw(5)_fbands(30)_win(NaN)ms_spec(gammatone)';
 

fn.load.file = sprintf(fn.load.file_template, data_type);
fn.load.fullfile = fullfile( fn.load.path, fn.load.file );
data = load(fn.load.fullfile);

spec_st = data.spec_st;
stim_st = data.stim_st;
tbl_data  = data.(sprintf('tbl_%s', data_type));
n_units   = height(tbl_data);     % total available units
duration_sec = 36;      % (sec) 
assert(duration_sec == spec_st.duration_ms * 1e-3,...
    '--> ERROR: You are using the wrong stimulus duration!');

aux.vprint(verbose, '--> [main_loopover_units.m] Loading file:\n\t...<%s>\n', fn.load.file);
aux.vprint(verbose, '-> data_type: %s\n', data_type);



%%
drr   = get_DRR_list_and_indices;
n_drr = drr.n_drr;                  % # DRRs of used 



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

   


%% Loop over all n_units for SELECTED DRR condition
bfcc_st = {};
tbl_BFcc = table();

for ix = 1:drr.n_drr
    %ix     = 1
    drr_ix = drr.ordered( ix );
    Sdrr   = spec_st.Sft{drr_ix};     % fry spectrogram
    n_freq = size(Sdrr, 1); 

    % Best-frequency correlation coefficient
    BF = nan(n_units, 1);
    CC  = nan(n_units, 1);
    P   = nan(n_units, 1);
    N = size(Sdrr,2);

    for k = 1:n_units     
        [BF(k), CC(k), P(k)] = BFcc(data.H(:,drr_ix,k), Sdrr, spec_st.f, false);      
    end
    
    % Add BF as columns into the measurement table
    % Create a new table with all the information of the measurements and save it
    %neuron = tbl_data.neuron;
    tbl_BFcc = [tbl_BFcc,...
        table(BF, 'VariableNames', {sprintf('BF%d',ix)}), ...
        table(CC, 'VariableNames', {sprintf('CC%d',ix)}), ...
        table(P, 'VariableNames', {sprintf('P%d',ix)}) ...
    ];     
end

neuron = tbl_data.neuron;
tbl_BFcc = [table(neuron), tbl_BFcc];

tbl_BFcc(1:5,:)

if save_results
    save([fn.load.fullfile, '_BFcc'], 'bfcc_st', 'tbl_BFcc', 'spec_st', 'stim_st');
end



%%
figure(11);
clf;
nbins = 25;
for k = 1:drr.n_drr
    drr_ix = drr.n_drr-k+1;
    subplot(drr.n_drr, 1, drr_ix);
    h = histogram(tbl_BFcc.(sprintf('CC%d',k)), nbins);
    h.FaceColor = aux.rpalette(k);
    ylabel(drr.labels{drr.ordered(k)});
end
title('Best Envelope Correlation');





