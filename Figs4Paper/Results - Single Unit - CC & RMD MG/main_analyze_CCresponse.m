% main_analyze_CCresponse.m
%
%

clc
fignum = 11;
verbose = 1;

setup_environment('../../');

duration_sec = 36;      % (sec) 

drr   = get_DRR_list_and_indices;
n_drr = drr.n_drr;                  % # DRRs of used 
idx_dry = drr.ordered(1);



%% Load data
%   Run [main_aggregate_MUA_data.m] again to update this file if needed
% 
fpath = load.path_to_data('_data');
fn_template = 'data_%s_(08-Jan-2021)_bw(5)_fbands(30)_win(NaN)ms_spec(gammatone).mat';
data_types = {'SU', 'MUA'};

clear data scores



%% SUs
% dtype           = data_types{1};
% data.su.fn      = fullfile(fpath, sprintf(fn_template, dtype));
% data.su.data    = load(su.fn);
% data.su.tbl     = su.data.(sprintf('tbl_%s', dtype));
% data.su.n_units = height(tbl_data);     % total available units
%     assert(duration_sec == spec_st.duration_ms * 1e-3,...
%         '--> ERROR: You are using the wrong stimulus duration!');

for q = 1:length(data_types)    
    dtype               = data_types{q};
    data.(dtype).fn     = fullfile(fpath, sprintf(fn_template, dtype));
    loaded_data         = load(data.(dtype).fn);
    data.(dtype).tbl    = loaded_data.(sprintf('tbl_%s', dtype));
    data.(dtype).n_units= height(data.(dtype).tbl);     % total available units
    assert(duration_sec == loaded_data.spec_st.duration_ms * 1e-3,...
        '--> ERROR: You are using the wrong stimulus duration!');
    
    H = loaded_data.H;
    spec_st = loaded_data.spec_st;
    n_units = data.(dtype).n_units;
    scores.(dtype).CC = nan(data.(dtype).n_units, n_drr);       % % CC(response dry & response)

    [sorted_list, tbl_BFcc] = find_best_unit_set('CC', 'fn', data.(dtype).fn); 
    
    for n = 1:n_units        
        % Find best frequency match
        [~, idx_bf] = min(abs(spec_st.f - tbl_BFcc.BF(n)));

        % Signal Envelope
        yenv = spec_st.Sft{ drr.dry }(idx_bf,:); 

        for k = 1:n_drr
            rv = drr.ordered(k);
            hn = H(:,rv,n);
            Rn = corrcoef(yenv, hn);
            scores.(dtype).CC(n,k) = Rn(1,2);
        end
    end
    
    % Rearrange in the desired order:
    scores.(dtype).CC = scores.(dtype).CC(sorted_list,:);
    
end



%% Correlation Coefficients between STIMULI
CCs = CC_stimuli( spec_st );



%% === SORTed ===
% CCer_SORTed: correlation coefficients between bandpass-envelope and response
figure(0+fignum);
clf;
linewidth = 5;
markersize = 25;
fontsize = 28;
fontsize_big = 32;
fontsize_bigger = 45;
ax = [];
h = [];

n2plot = 103;

for q = 1:length({'SU', 'MUA'})
    ax(q) = subplot(1,2,q);
    dtype = data_types{q};
    CCk = scores.(dtype).CC(1:n2plot,:);
    h = aux.violinplot(CCk, drr.labels(drr.ordered),...
        'ShowMean', true,...
        'ViolinAlpha', 0.5);    
        %, 'ShowData', false);
    for jj = 1:n_drr
        h(jj).ViolinColor = aux.rpalette(jj);
    end
    set(ax, 'FontSize', fontsize);
    hold on
    plth = plot(CCs, 'sk:', 'MarkerSize', markersize, 'MarkerFaceColor', 'k');
    hold off
    title(ax(q), sprintf('%d %s', size(CCk,1), dtype), 'FontSize', fontsize_bigger);
    xlabel(ax(q), 'DRR', 'FontSize', fontsize_big);

end

ylabel(ax(1), 'CC (Envelope-to-Response)', 'FontSize', fontsize_big);
linkaxes(ax);
ylim(ax(1), [-0.15, 1.1]);
set(ax(2), 'YTickLabel', '');
xlim([0.6, 5.4]);
set(ax(1), 'Position', [0.0970, 0.1100, 0.3928, 0.8038]);
set(ax(2), 'Position', [0.5279, 0.1100, 0.3928, 0.8038]);








