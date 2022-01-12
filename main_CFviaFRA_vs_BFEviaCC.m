%
% main_FRA.m.m
% 
% Description:
%   Reads all FRAs from the raw data.
% 
%

% viewer.plot_FRA(S);


clc
% clear all

fignum = 11;
verbose     = 0;

setup_environment('../');
data_path = load.path_to_data('raw');
assert(7 == exist(data_path, 'dir'), sprintf('DirectoryError: no such directory: %s !', data_path));

warning('off');


%% Load data
%   Run [main_aggregate_MUA_data.m] again to update this file if needed
% 
data_type   = 'SU';       % {'SU', MUA'}
fn.load.path= load.path_to_data('_data');
data_type   = upper(data_type);

fn.load.file_template = 'data_%s_(08-Jan-2021)_bw(5)_fbands(30)_win(NaN)ms_spec(gammatone)';
% fn.load.file_template = 'data_%s_(05-Jan-2022)_bw(5)_fbands(30)_win(NaN)ms_spec(gammatone-SYNC)';
% fn.load.file_template = 'data_%s_(10-Jan-2022)_bw(5)_fbands(30)_win(NaN)ms_spec(gammatone-SYNC)';

fn.load.file = sprintf(fn.load.file_template, data_type);
fn.load.fullfile = fullfile( fn.load.path, [fn.load.file, '.mat'] );
dummy = load(fn.load.fullfile, 'tbl_SU');
tbl_SU = dummy.tbl_SU;
n_units = height(tbl_SU);



%% Get the BFcc:
% dir([fn.load.fullfile(1:end-4), '*'])
bfcc_filename = fullfile( fn.load.path, [fn.load.file, '_BFcc.mat'] );
dummy = load(bfcc_filename);
tbl_BFcc = dummy.tbl_BFcc;




%% 
syncchan = 0;   % Impale notation for syncronized break/new repetition.   
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
        
        
        % Get a list of all channel numbers/spike numbers
        all_spikechan = spiketools.get_all_spikechan(S, syncchan);
        if iscell(all_spikechan), all_spikechan = all_spikechan{:}; end

        if 0 == sum( all_spikechan == spikechan )
            spikechan = 1;
        end
        
        % Calc FRA for each of the available spikes
        fra = arrayfun(@(SPK) medit.FRA(S, SPK, syncchan), spikechan, 'UniformOutput', false) ;
        assert(1==length(fra), sprintf('length(fra) > 1 !'));
        cf_m = fra{1}.CF;
        
        fprintf('\n');
        fprintf(' (%d) spikechan: %.2f\n', k, spikechan);
        fprintf(' (%d) cf       : %.2f Hz\n', k, cf_m);

        
        if false & ~isnan(cf_m)    
            figure(1)
            hold on
            hk = plot(1e-3*cf_m, 1e-3*tbl_BFcc.BF(k), 'k+', 'MarkerSize', 20, 'LineWidth', 2);
            ax = viewer.plot_FRA(S, 'spikechan', spikechan, 'figh', figure(2));
            %pause(1.0);            
            set(ax, 'XTickLabel', '', 'YTickLabel', '');
            xlabel('');
            ylabel('');
            %caxis([0, 80]);
            title('');
            
            delete(hk)
        end 
        
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




%% Get the latest CF
cf = nan(n_units,1);

for k = 1:n_units
    cf_krow = tbl_fra(k,2:end).Variables;

    % Option #1
    %ix = find(~isnan(cf_krow), 1, 'last');
    
    % Option #2
    bfcc = tbl_BFcc.BF(k);
    [~, ix] = min(abs(cf_krow - bfcc).^2);
    
    if all(isnan(cf_krow))
        continue;
    end
    cf(k) = cf_krow(ix);
end




%% Remove "bad" FRA readings
ix_bad_fra = [3, 18, 23:24, 30, 32, 37:38, 40, 42, 45:46, 48:50, 53:56, 58, 60, 68:71, ...
    73:76, 79, 81:83, 85:89, 91, 93, 94, 97:98, 101:102];

cf(ix_bad_fra) = nan;
tbl_fra.cf_final = cf;
tbl_fra.cf_final(51) = tbl_fra.cf4(51);

% [tbl_fra(~isnan(tbl_fra.cf_final),:), tbl_BFcc(~isnan(tbl_fra.cf_final),'BF')]



%%
un  = 1e-3;

ix  = ~isnan(tbl_fra.cf_final);
% ix  = 1 ==( ~isnan(tbl_fra.cf_final) .* (tbl_fra.cf_final < 8000) );

cf  = un * tbl_fra.cf_final( ix );
bfcc= un * tbl_BFcc.BF( ix );
R   = 500*abs(tbl_BFcc.R( ix ));


%% CFs under 8k Hz
ix_ = ix & (tbl_fra.cf_final < 8000);
cf_ = un * tbl_fra.cf_final( ix_ );
bfcc_= un * tbl_BFcc.BF( ix_ );

rr = corrcoef(cf_, bfcc_);
rr = rr(1,2);
fprintf('\n- corrcoef(CFs, BFcc): %.2f (R2)\n', rr);
fprintf('- # of CFs: %d \n', sum(ix));
fprintf('- # of CFs: %d (regression)\n', sum(ix_));



%% R square with corrcoef
fontsize = 24;
markersize = 52;

figure(5);
clf;

% Add a bit of jitterness
x_rnd = 0.05*randn(sum(ix), 1);
y_rnd = 0.05*randn(sum(ix), 1);

%plot(1e-3*tbl_fra.cf_final, 1e-3*tbl_BFcc.BF, '.', 'MarkerSize', markersize);
%h = scatter(1e-3*tbl_fra.cf_final, 1e-3*tbl_BFcc.BF, 200*abs(tbl_BFcc.R), 'filled');
h = scatter(cf + x_rnd, bfcc + y_rnd, R,...
    'filled', 'MarkerFaceAlpha',.65 );

xylabels = xlabel('CF Frequency (Hz)');
xylabels(2) = ylabel('BFcc Frequency (Hz)');
hold on
plot(un*[0,18e3], un*[0,18e3], ':');
hold off
xlim(un*[0, 16e3]);
ylim(un*[0, 9e3]);
axis square
aux.vline(un*8000, 'LineWidth', 1);
aux.hline(un*8000, 'LineWidth', 1);



% Fit linear regression line with OLS
b = [ones(size(cf_,1),1) cf_]\bfcc_;
% Use estimated slope and intercept to create regression line
rline = [ones(size(cf_,1),1) cf_]*b;    % regression line
hold on
hreg = plot(sort(cf_), sort(rline), '--k', 'LineWidth', 4);
hold off

hreg.DisplayName = sprintf('Regression Line $(y = %0.2f\\cdot x + %0.2f)$',b(2),b(1));
legend(hreg);

title( sprintf('%d Units $(r^2: %.2f)$', sum(~isnan(cf)), rr), 'fontsize', fix(1.8*fontsize ));


set(gca, 'fontsize', fix(1.4*fontsize));
set(xylabels, 'fontsize', fix(1.4*fontsize));





%% R square (another way to compute)
% Add a regression line
%
% See: https://www.mathworks.com/matlabcentral/answers/478999-how-to-show-r-square-correlation-and-rmse-on-a-scatterplot
%
y = bfcc_;

% RMSE between regression line and y
RMSE = sqrt(mean((y-rline).^2));
% R2 between regression line and y
SS_X = sum((rline-mean(rline)).^2);
SS_Y = sum((y-mean(y)).^2);
SS_XY = sum((rline-mean(rline)).*(y-mean(y)));
R_squared = SS_XY/sqrt(SS_X*SS_Y);
fprintf('RMSE: %0.2f | R2: %0.2f\n',RMSE,R_squared)




%% log-scale
%{
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
ax = gca;
ax.FontSize = fontsize;
set(xylabels, 'fontsize', fix(1.4*fontsize));
title(sprintf('%d Units', sum(~isnan(cf))), 'fontsize', fix(1.8*fontsize))
%}



