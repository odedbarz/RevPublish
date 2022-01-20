% main_neural_noise.m
%
%
%


clc
fignum = 11;
verbose = 1;

setup_environment('..');

duration_sec = 36;      % (sec) 

drr   = get_DRR_list_and_indices;
n_drr = drr.n_drr;                  % # DRRs of used 
idx_dry = drr.ordered(1);



%% Load data
%   Run [main_aggregate_MUA_data.m] again to update this file if needed
% 
fn.load.path = load.path_to_data('_data');
fn.load.file_template = 'data_%s_(13-Jan-2022)_bw(5)_fbands(30)_win(NaN)ms_spec(gammatone-only)';
% data_types = {'SU', 'MUA'};
data_type = 'MUA';

file2load = fullfile(fn.load.path, [sprintf(fn.load.file_template, data_type), '_ALL-MEAS.mat']);

dummy = load(file2load);

spec_st     = dummy.spec_st;
stim_st     = dummy.stim_st;
T = dummy.tbl_MUA_all;    % tbl_MUA_all
clear dummy



%% Calc coefficient of variation (CV) 
n_units = height(T);
dim = 1;
CV = nan(n_units, drr.n_drr);

for k = 1:n_units
    % Reorder the (sub)table of measurements according the DRRs
    Tk = T.all_meas{k};
    if height(Tk) < 5
        n_units = n_units -1;
        continue;
    end
        
    Tk = Tk(drr.ordered,:);
    
    %clf
    
    for dr = 1:drr.n_drr
        mu = mean( Tk.single_meas{dr}, dim );
        sig= std( Tk.single_meas{dr}, [], dim ); 
        CV(k, dr) = mean( sig./mu );
        
        %plot( hist( Tk.single_meas{dr}(:), 100, 'r' ) );
        %hold on
        %legend(drr.labels(drr.ordered))
    end
    
end



%%
figure(fignum+1);
clf;
boxplot(CV, drr.labels(drr.ordered));
ylabel('$<\sigma/\mu>$');
xlabel('DRR');
set(gca, 'FontSize', 24);
title( sprintf('%s (%d units)', data_type, n_units) );



