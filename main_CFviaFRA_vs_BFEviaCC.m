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
verbose= 0;

setup_environment('../');
data_path = load.path_to_data('raw');
assert(7 == exist(data_path, 'dir'), sprintf('DirectoryError: no such directory: %s !', data_path));

warning('off');


%% Load data
%   Run [main_aggregate_MUA_data.m] again to update this file if needed
% 
data_type   = 'MUA';       % {'SU', MUA'}
fn.load.path= load.path_to_data('_data');
data_type   = upper(data_type);

fn.load.file_template = 'data_%s_(13-Jan-2022)_bw(5)_fbands(30)_win(NaN)ms_spec(gammatone-only)';


fn.load.file = sprintf(fn.load.file_template, data_type);
fn.load.fullfile = fullfile( fn.load.path, [fn.load.file, '.mat'] );
tbl_name = sprintf('tbl_%s',data_type);
dummy = load(fn.load.fullfile, tbl_name);

T = dummy.(tbl_name);
% tbl_data  = data.(sprintf('tbl_%s', data_type));

n_units = height(T);



%% Get the BFcc:
% dir([fn.load.fullfile(1:end-4), '*'])
bfcc_filename = fullfile( fn.load.path, [fn.load.file, '_BFcc.mat'] );
dummy = load(bfcc_filename);
tbl_BFcc = dummy.tbl_BFcc;



%% Order the units
sort_type = 'SPK';  % {'RND', 'SVD', 'FILE', 'SPK', 'NOSPK'}
ix_spk = find_best_unit_set(sort_type,...
    'fpath', load.path_to_data('_data'),...
    'fn_template', fn.load.file_template, ...
    'data_type', data_type);
spk_list = find(ix_spk);



%% 
syncchan = 0;   % Impale notation for syncronized break/new repetition.   
max_repeated_trials = 4;   % max num of channels in this data (of LMAs)
var_names = arrayfun(@(x) sprintf('cf%d', x), 1:max_repeated_trials, 'UniformOutput', false);
tbl_fra = array2table(nan(n_units, max_repeated_trials), 'VariableNames', var_names);    


for k = 1:n_units
    if ~any(k == spk_list)
        continue;
    end
    session = T.expName{k};
    session_path = fullfile(data_path, session);
    assert(7 == exist(session_path, 'dir'), sprintf('DirectoryError: no such directory: %s !', data_path));
    
    files = dir( session_path );
    files = {files.name};
    TF1 = contains(files, sprintf('%s-%d-', session, T.unit(k)));
    TF2 = contains(files, 'FRA');
    fra_idx = find( TF1 & TF2 );
    if isempty(fra_idx)
        fra_idx = find( TF2 );  % in these sessions I did only one FRA test!
    end        
    num_fra_sessions = length(fra_idx);  % number of FRA files for a given session/directory
    spikechan = T.spikechan(k);        

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
        %fra = arrayfun(@(SPK) medit.FRA(S, SPK, syncchan), spikechan, 'UniformOutput', false) ;
        %assert(1==length(fra), sprintf('length(fra) > 1 !'));
        fra =  medit.FRA(S, spikechan, syncchan);
%         CF = medit.calc_CF(S, fra.rates, 10, 1, 5*2.5);
        cf_m = fra.CF;
        
        fprintf('\n');
        fprintf(' (%d) spikechan: %d\n', k, spikechan);
        fprintf(' (%d) cf       : %.2f Hz\n', k, cf_m);
        
%                 viewer.plot_FRA(S, 'spikechan', spikechan);
%                 fprintf(' - spikes: %d\n', sum(sum(fra.spike_count)));

        if sum(fra.spike_count(:)) < 20, continue; end
        
        if isempty(cf_m), continue; end
        
        
        tbl_fra{k, c} = cf_m;
        c = c + 1;
    end
    
    
    if verbose
        fprintf('\t(%d) (neuron: %d; spikechan: %d)\n', k, T.neuron(k), spikechan);
        disp( T(1:k,:) );
        disp( tbl_fra(1:k,:) );
    end
end

% Add the neuron's numbers too
tbl_fra = [T(:,'neuron'), tbl_fra];




%% Get the latest CF
cf = nan(n_units,1);

for k = 1:n_units
    cf_k_row = tbl_fra(k,2:end).Variables;

    % Option #1
    %ix = find(~isnan(cf_k_row), 1, 'last');
    
    % Option #2
    bf1 = tbl_BFcc.BF1(k);
    [~, ix] = min(abs(cf_k_row - bf1).^2);
    
    if all(isnan(cf_k_row))
        continue;
    end
    cf(k) = cf_k_row(ix);
end






%% Manually removing "bad" FRA that the algorithm didn't catch
% ix_bad_fra = [3, 18, 23:24, 30, 32, 37:38, 40, 42, 45:46, 48:50, 53:56, 58, 60, 68:71, ...
%     73:76, 79, 81:83, 85:89, 91, 93, 94, 97:98, 101:102];
ix_bad_fra_su = [3, 23, 43, 45:46, 51, 53:57, 60, 70, 74, 76, 81:94]; %, 101:102];
ix_bad = spk_list(ix_bad_fra_su);

cf(ix_bad) = nan;
tbl_fra.cf_final = cf;
tbl_fra.cf_final(51) = tbl_fra.cf4(51);

% [tbl_fra(~isnan(tbl_fra.cf_final),:), tbl_BFcc(~isnan(tbl_fra.cf_final),'BF')]



% tbl_fra.cf_final = cf;  '@@@@@@@@@@@ DEBUG @@@@@@@@@@@'


%%
un  = 1e-3;

% ix  = ~isnan(tbl_fra.cf_final);
ix  = ~isnan(tbl_fra.cf_final); % & ...
%     tbl_BFcc.CC1 > median(tbl_BFcc.CC1) & ...
%     (tbl_fra.cf_final < 8000);
% ix  = 1 ==( ~isnan(tbl_fra.cf_final) .* (tbl_fra.cf_final < 8000) );

cf  = un * tbl_fra.cf_final( ix );
bf1= un * tbl_BFcc.BF1( ix );
cc1   = 500*abs(tbl_BFcc.CC1( ix ));
% cc1 = 200;


%% CFs under 8k Hz
ix_ = ix & (tbl_fra.cf_final < 8000);

% % weights:
% w = 1/sum(ix_)*ones(sum(ix_), 1);    
% % w = abs(tbl_BFcc.CC1(ix_));  % CC weights
% 
cf_ = un * tbl_fra.cf_final( ix_ );
% cf_mean = cf_(:) .* w/sum(w);
% 
bfcc_= un * tbl_BFcc.BF1( ix_ );
% bfcc_mean = bfcc_(:) .* w/sum(w);
% 
% wcovxy = (w'*((cf_ - cf_mean) .* (bfcc_ - bfcc_mean)))/sum(w);
% wcovx = (w'*(cf_ - cf_mean).^2)/sum(w);
% wcovy = (w'*(bfcc_ - bfcc_mean).^2)/sum(w);
% cc = wcovxy/sqrt( wcovx * wcovy );

[cc, pv] = corrcoef(cf_, bfcc_);
cc = cc(1,2);
pv = pv(1,2);

fprintf('\n- corrcoef(CFs, BFcc): %.2f (R2), p: %g\n', cc, pv);
fprintf('- # of units: %d (total)\n', sum(ix));
fprintf('- # of units: %d (for regression)\n', sum(ix_));



%% R square with corrcoef
fontsize = 24;
markersize = 52;

figure(10 + strcmpi('MUA', data_type)*5);
clf;

% Add a bit of jitterness
x_rnd = 0.05*randn(sum(ix), 1);
y_rnd = 0.05*randn(sum(ix), 1);

h = scatter(cf + x_rnd, bf1 + y_rnd, cc1, ...
    'filled', 'MarkerFaceAlpha',.65, 'DisplayName', data_type );

if strcmpi('SU', data_type)
    h.MarkerFaceColor = aux.rpalette(1);
else
    h.MarkerFaceColor = aux.rpalette(2);
end

xylabels = xlabel('CF Frequency (Hz)');
xylabels(2) = ylabel('BFcc Frequency (Hz)');
hold on
plot(un*[0,18e3], un*[0,18e3], ':');
hold off
xlim(un*[0, 16e3]);
ylim(un*[0, 9e3]); 
% axis square
aux.vline(un*8000, 'LineWidth', 1);
aux.hline(un*8000, 'LineWidth', 1);


title( sprintf('%s (%d Units, $r^2: %.2f$)', data_type, sum(~isnan(cf)), cc), 'fontsize', fix(1.8*fontsize ));
set(gca, 'fontsize', fix(1.4*fontsize));
set(xylabels, 'fontsize', fix(1.4*fontsize));




%% Fit linear regression line with OLS
b = [ones(size(cf_,1),1) cf_]\bfcc_;
% Use estimated slope and intercept to create regression line
rline = [ones(size(cf_,1),1) cf_]*b;    % regression line
hold on
hreg = plot(sort(cf_), sort(rline), '--k', 'LineWidth', 4);
hold off

hreg.DisplayName = sprintf('Regression Line $(y = %0.2f\\cdot x + %0.2f)$',b(2),b(1));
legend(hreg);








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
fprintf('%s:\n - RMSE: %0.2f\n - R2  : %0.2f\n',data_type, RMSE,R_squared)




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



