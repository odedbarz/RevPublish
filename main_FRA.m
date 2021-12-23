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

% setup_environment('../');
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
        cf = fra{1}.CF;
        if isempty(cf), continue; end
        tbl_fra{k, c} = cf;
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



%%
figure(1);
clf;
plot(tbl_fra.cf1, tbl_BFcc.BF, '.', 'MarkerSize', 32);
xlabel('CF Frequency (Hz)');
ylabel('BFcc Frequency (Hz)');
hold on
for k = 2:max_repeated_trials
    plot(tbl_fra.(sprintf('cf%d',k)), tbl_BFcc.BF, 'x');
end
plot([0,18e3], [0,18e3], ':');
hold off
max_freq = 14e3;
xlim([0, max_freq]);
ylim([0, max_freq]);
axis square
aux.vline(8000, 'LineWidth', 1);
aux.hline(8000, 'LineWidth', 1);







