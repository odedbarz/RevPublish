%
% main_raw_MUA_stats.m
%
% Description:
% Check for statistics within raw measurements
%
%
%

clc
fignum = 11;
verbose = 1;

setup_environment;


%% Paths
% Path to the Impale's data
path_root_mat = load.path_to_data('Impale_data');

% Path to the RAW data
path_root_raw = load.path_to_data('raw');



%% Load the measurement's table
if ~exist('tbl_main', 'var')
    % Save time, load this table only once 
    tbl_impale = readtable([path_root_mat, 'LabNoteBook_Reverb.xlsx'], 'Sheet', 'Spch');
    
    % Row is a duplicate variable name used by MATLAB
    tbl_impale.Properties.VariableNames{'Row'} = 'neuron';
end
    


%% *** Get the MUA data into one big matrix ***
% ALways load stimuli with the same duration
'Loads stimuli with the SAME duration'
duration_sec= 36;   % (sec) stimulus duration to use 
duration_ms = units.sec2ms( duration_sec );
binwidth    = 1         % (ms) binwidth of the resulted spectrogram 
fs          = 1/(1e-3*binwidth);         % (Hz)

if verbose
    aux.cprintf('UnterminatedStrings', '\n    Stimulus Parameters:\n');    
    aux.cprintf('UnterminatedStrings', '--> Selected DURATION: %g (sec)\n', duration_sec);    
    aux.cprintf('UnterminatedStrings', '--> binwidth         : %g ms\n', binwidth);
end

% A 3D matrix to hold all the MUA responses
drr     = get_DRR_list_and_indices;
n_drr 	= drr.n_drr;  
fs_raw  = 10e3;     % (Hz) that's the frequency that I save in the raw table files
fs_mua  = fs;       % (Hz) that's the new downsampled frequency of the MUA data
n_smp   = fs_mua * duration_sec;


% Get all rows (measurements) with desired duration
neuron_list = find(tbl_impale.duration_sec == duration_sec);

% Get all units (rows) with duration_sec
tbl_MUA = tbl_impale(neuron_list, :);
n_rows  = height(tbl_MUA);
H_labels= zeros(n_rows, n_drr);

% Initialize the statistical structure
clear stats
stats.drr_labels = drr.labels(drr.ordered);
stats.CC.CC  = cell(n_drr,1);
stats.CC.grp = cell(n_drr,1);
stats.CC.n   = zeros(n_drr, 1);
stats.CC.mean = zeros(n_drr, 1);
stats.CC.median = zeros(n_drr, 1);
stats.CC.std = zeros(n_drr, 1);
stats.CC.mad = zeros(n_drr, 1);            
stats.snr = zeros(n_drr, 1);



%% Prepare the MUA measurements to match the spectrograms
aux.vprint(verbose, '\n-> Starting main loop...\n');
for k = 1:n_rows
    tbl_neuron_k = tbl_MUA(k,:);
    S = load.response( tbl_neuron_k, duration_sec );
    if isempty(S)
        aux.vprint(verbose, '--> Skipping measurement (%d/%d)...\n', k, n_rows);
        continue;
    elseif 0==mod(k-1,10)
        aux.vprint(verbose, '--> Loading Imaple structure (%d/%d)\n', k, n_rows);
    end
    S  = S{1};     % cell --> struct

    % Loads the raw data
    fn.raw = fullfile(path_root_raw, S.info.rawDataDir, S.info.rawFileName);
    assert(~isempty(dir(fn.raw)), '--> ERROR: can''t find this file!!\n');
    raw_st = load(fn.raw, '-mat');
    assert(fs_raw == raw_st.sr);

    % Calculate MUA 
    [~, ~,  tbl_kth_meas] = calc_raw_means(raw_st.tbl, raw_st.sr, binwidth, duration_sec, tbl_neuron_k.spikechan, 'MUA');    
        
    % Skip uncomplete measurements
    if n_drr >= height(tbl_kth_meas)
        continue;
    end
        
    % Make sure to use ALL measurements (including the second DRY
    % measurement)
    %assert( isempty(tbl_kth_meas.single_meas{end}),...
    %    '--> Last measurement (DRY) is not empty -- use it!');
    if ~isempty(tbl_kth_meas.single_meas{end})
        Hdry = tbl_kth_meas.single_meas{drr.ordered(1)};
        assert(1.5 == tbl_kth_meas.Dist(drr.ordered(1)));
        assert(100 == tbl_kth_meas.Reverb(drr.ordered(1)));
        
        Hdry2 = tbl_kth_meas.single_meas{drr.sortby(end)};
        assert(3.0 == tbl_kth_meas.Dist(drr.sortby(end)));
        assert(100 == tbl_kth_meas.Reverb(drr.sortby(end)));
        
        tbl_kth_meas.single_meas{drr.ordered(1)} = [Hdry, Hdry2];
    end
        
    % Get all measurements with N_DRR (or more); don't add incomplete measurements 
    H_labels(k, :) = arrayfun(@(N) ~isempty(tbl_kth_meas.single_meas{N}), 1:n_drr);    
    if n_drr > sum(H_labels(k,1:n_drr),2)
        continue;
    end
    
    % Make sure that the order of the DRR conditions is always the "right" one
    assert( all(drr.revb(1:n_drr) == tbl_kth_meas.Reverb(1:n_drr)'),...
        '!ERROR: Wrong order of measurements!' )    
    
    % Calculate pairwise-distance between all vectors. 
    for dr = 1:n_drr
        drr_idx = drr.ordered(dr);    
        Y = tbl_kth_meas.single_meas{drr_idx};
        
        % ### Correlation between trials ###
        Mcc = corr(Y, 'Type', 'Pearson');
        
        % -> Get the upper off-diagonals values. But distance matrix Z must 
        %    be square with 0 along the diagonal.
        CCk = 1.0 - squareform(1 - Mcc, 'tovector');        
        
        % Add to the right row
        stats.CC.n(dr)   = 1 + stats.CC.n(dr); 
        len_CC        = length(CCk);
        stats.CC.grp{dr} = [stats.CC.grp{dr}, stats.CC.n(dr)*ones(1, len_CC)];
        stats.CC.CC{dr}  = [stats.CC.CC{dr}, CCk];
        
        stats.CC.median(dr, stats.CC.n(dr)) = median(CCk);
        stats.CC.mean(dr, stats.CC.n(dr))   = mean(CCk);
        stats.CC.std(dr, stats.CC.n(dr))    = std(CCk);
        stats.CC.mad(dr, stats.CC.n(dr))    = mad(CCk);
        
        
        % ### Signal-to-Noise between trials ###
        Yavg = mean(Y, 2);
        stats.snr(dr, stats.CC.n(dr)) = median( var(Yavg)./var(Yavg-Y+eps) );
        
    end
    tbl_MUA.drr{k} = CCk;
       
end

arrayfun(@(N) assert( length(stats.CC.CC{N}) == length(stats.CC.grp{N}),...
    'ERROR: something is wrong! number of GROUPs and number of CCs must be the same!'),...
    1:n_drr);



%% Get all rows (measurements) that have full session (all DRR conditions)
valid_neuron_idx = n_drr == sum(H_labels(:,1:n_drr),2);
tbl_MUA = tbl_MUA(valid_neuron_idx, :);



%% Save the data
% %{
fprintf('SAVE the analysis!');

fn.path = load.path_to_data('Stats');
fn.file = sprintf('stats_raw(MUA)_(%s)_BW(%g)ms_duration(%d)sec_units(%d)', ...
    date, binwidth, duration_sec, height(tbl_MUA));
fn.save = fullfile(fn.path, fn.file);

fprintf('\nSaving data at:\n');
disp(fn)
save(fn.save, 'stats', 'tbl_MUA');
%}




