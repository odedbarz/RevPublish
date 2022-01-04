%
% main_FRA.m.m
% 
% Description:
%   Reads all FRAs from the raw data.
% 
%

% viewer.plot_FRA(S);


clc
clear all

fignum = 11;
plot_things = 0;
verbose     = 0;

setup_environment('../');
data_path = load.path_to_data('raw');
assert(7 == exist(data_path, 'dir'), sprintf('DirectoryError: no such directory: %s !', data_path));



%% Load data
%   Run [main_aggregate_MUA_data.m] again to update this file if needed
% 
data_type   = 'SU';       % {'SU', MUA'}
fn.load.path= load.path_to_data('_data');
data_type   = upper(data_type);

fn.load.file_template = 'data_%s_(08-Jan-2021)_bw(5)_fbands(30)_win(NaN)ms_spec(gammatone)';
% fn.load.file_template = 'data_%s_(01-Nov-2021)_bw(1)_fbands(30)_win(NaN)ms_spec(gammatone)';
% fn.load.file_template = 'data_%s_(08-Nov-2021)_bw(5)_fbands(30)_win(NaN)ms_spec(gammatone)';
% fn.load.file_template = 'data_%s_(10-Dec-2021)_bw(100)_fbands(30)_win(NaN)ms_spec(gammatone)'

fn.load.file = sprintf(fn.load.file_template, data_type);
fn.load.fullfile = fullfile( fn.load.path, [fn.load.file, '.mat'] );
dummy = load(fn.load.fullfile, 'tbl_SU');
tbl_SU = dummy.tbl_SU;
n_units = height(tbl_SU);



%% 
syncchan = 0;   % Impale notation for syncronized break/new repetition.   



%%
max_repeated_trials = 4;   % max num of channels in this data (of LMAs)
var_names = arrayfun(@(x) sprintf('cf%d', x), 1:max_repeated_trials, 'UniformOutput', false);
tbl_fra = array2table(nan(n_units, max_repeated_trials), 'VariableNames', var_names);    

for k = 1:n_units
    session = tbl_SU.expName{k};
    session_path = fullfile(data_path, session);
    assert(7 == exist(session_path, 'dir'), sprintf('DirectoryError: no such directory: %s !', data_path));
    
    files = dir( session_path );
    files = {files.name};
    TF1 = contains(files, sprintf('%s-%d-', session, tbl_SU.unit(k)));
    TF2 = contains(files, 'FRA');
    fra_idx = find( TF1 & TF2 );
    if isempty(fra_idx)
        fra_idx = find( TF2 );  % in these sessions I did only one FRA test!
    end        
    num_fra_sessions = length(fra_idx);  % number of FRA files for a given session/directory
    spikechan = tbl_SU.spikechan(k);        

    % In case of multiplications, saves the LAST available FRA 
    c = 1;  % counter for the table
    for m = 1:num_fra_sessions
        fn_fra = fullfile(session_path, files{fra_idx(m)});
        assert(2 == exist(fn_fra, 'file'), sprintf('FilenameError: no such file: %s !', fn_fra));
        S = load(fn_fra);
        if plot_things
            viewer.plot_FRA(S, fignum, syncchan);
        end
        
        % Get a list of all channel numbers/spike numbers
        all_spikechan = spiketools.get_all_spikechan(S, syncchan);
        if iscell(all_spikechan), all_spikechan = all_spikechan{:}; end

        if 0 == sum( all_spikechan == spikechan )
            spikechan = 1;
        end
        
        % Calc FRA for each of the available spikes
        fra = arrayfun(@(SPK) medit.FRA(S, SPK, syncchan), spikechan, 'UniformOutput', false) ;
        cf_m = fra{1}.CF;
        
        fprintf('\n');
        fprintf(' (%d) spikechan: %.2f Hz\n', k, spikechan);
        fprintf(' (%d) cf       : %.2f Hz\n', k, cf_m);

%         viewer.plot_FRA(S, 'spikechan', spikechan);

%         if 1000 > cf_m
% %             viewer.plot_FRA(S);
%             viewer.plot_FRA(S, 'spikechan', spikechan);
%             [];
%         end
        
        if isempty(cf_m), continue; end
        tbl_fra{k, c} = cf_m;
        c = c + 1;
    end
    
    
    if verbose
        fprintf('\t(%d) (neuron: %d; spikechan: %d)\n', k, tbl_SU.neuron(k), spikechan);
        disp( tbl_SU(1:k,:) );
        disp( tbl_fra(1:k,:) );
    end
end

% Add the neuron's numbers too
tbl_fra = [tbl_SU(:,'neuron'), tbl_fra];


%% Get the BFcc:
% dir([fn.load.fullfile(1:end-4), '*'])
bfcc_filename = fullfile( fn.load.path, [fn.load.file, '_BFcc.mat'] );
dummy = load(bfcc_filename);
tbl_BFcc = dummy.tbl_BFcc;



%% Get the latest CF
cf = nan(n_units,1);

for k = 1:n_units
    cf_krow = tbl_fra(k,2:end).Variables;
    ix = find(~isnan(cf_krow), 1, 'last');
    cf(k) = cf_krow(ix);
end




%% Remove "bad" FRA readings
cf([18, 23, 46, 49, 50, 56, 74, 94]) = nan;
tbl_fra.cf_final = cf;



%%
fontsize = 24;
markersize = 52;

figure(1);
clf;
plot(1e-3*tbl_fra.cf_final, 1e-3*tbl_BFcc.BF, '.', 'MarkerSize', markersize);
% loglog(tbl_fra.cf_final, tbl_BFcc.BF, '.', 'MarkerSize', 32);
xylabels = xlabel('CF Frequency (Hz)');
xylabels(2) = ylabel('BFcc Frequency (Hz)');
hold on
plot(1e-3*[0,18e3], 1e-3*[0,18e3], ':');
hold off
xlim(1e-3*[0, 16e3]);
ylim(1e-3*[0, 9e3]);
axis square
aux.vline(1e-3*8000, 'LineWidth', 1);
aux.hline(1e-3*8000, 'LineWidth', 1);
ax = gca;
ax.FontSize = fontsize;
set(xylabels, 'fontsize', fix(1.4*fontsize));
title(sprintf('%d Units', sum(~isnan(cf))), 'fontsize', fix(1.8*fontsize))


%% log-scale
figure(5);
clf;
loglog(tbl_fra.cf_final, tbl_BFcc.BF, '.', 'MarkerSize', markersize);
xylabels = xlabel('CF Frequency (Hz)');
xylabels(2) = ylabel('BFcc Frequency (Hz)');
hold on
plot([100,18e3], [100,18e3], ':');
hold off
% max_freq = 14e3;
xlim([1e3, 16e3]);
ylim([100, 9e3]);
axis square
% aux.vline(8000, 'LineWidth', 1);
% aux.hline(8000, 'LineWidth', 1);
ax = gca;
ax.FontSize = fontsize;
set(xylabels, 'fontsize', fix(1.4*fontsize));
title(sprintf('%d Units', sum(~isnan(cf))), 'fontsize', fix(1.8*fontsize))




