%
% plot_est_script.m
%


legend_FS = 13;     % legend's FontSize
legend_location = 'northeast';

% Whitening sphere
zca = @(x) (x-mean(x))./std(x);

t_test = spec_st.t(test_st.idx);    % (sec)
t_test = t_test - t_test(1);

clear ax
ax.mua = [];
ax.lfp = [];
ax.psth = [];
len_mua = length(mua_list);

ax.all = aux.tight_subplot(8, 3);
ax.mua = flipud( ax.all(1:3:end) );
ax.lfp = flipud( ax.all(2+(1:3:end)) );
ax.psth = flipud( ax.all(1+(1:3:end)) );

score.cc.mua = nan(1,len_mua);
score.cc.psth = nan(1,len_mua);
score.cc.lfp = nan(1,len_mua);

for k = 1:len_mua
    % MUA
    %ax.mua = [ax.mua, subplot(len_mua, 3, 3*len_mua-1 - (3*k-2))];  
    axes( ax.mua(k) );
    h = plot(t_test, zca([mua_list{k}.r_test, mua_list{k}.r_est]));
    axis(ax.mua(k), 'tight');
    CC = corrcoef(mua_list{k}.r_test, mua_list{k}.r_est);
    score.cc.mua(k) = CC(1,2);
    set(get(get(h(1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off')
    %legend_h = legend(sprintf('$\\hat{r}(t)$ (CC: %.3f)',  score.cc.mua(k)));
    legend_h = legend(sprintf('%.3f',  score.cc.mua(k)), 'Location', legend_location);
    legend_h.FontSize = legend_FS;
    
    % Spikes (PSTH)
    %ax.psth = [ax.psth, subplot(len_mua, 3, 3*len_mua - (3*k-2))];
    axes( ax.psth(k) );
    if ~isempty(psth_list{k})
        h = plot(t_test, zca([psth_list{k}.r_test, psth_list{k}.r_est]));
        axis(ax.psth(k), 'tight');
        CC = corrcoef(psth_list{k}.r_test, psth_list{k}.r_est);
        score.cc.psth(k) = CC(1,2);
        set(get(get(h(1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off')
        %legend_h = legend('$r(t)$', sprintf('$\\hat{r}(t)$ %.3f',  score.cc.psth(k)));
        legend_h = legend(sprintf('%.3f',  score.cc.psth(k)), 'Location', legend_location);
        legend_h.FontSize = legend_FS;        
    else
        set(ax.psth(k), 'Visible', 'Off');
    end
    
    
    % LFP
    %ax.lfp = [ax.lfp, subplot(len_mua, 3, 3*len_mua+1 - (3*k-2))];  
    axes( ax.lfp(k) );
    h = plot(t_test, zca([lfp_list{k}.r_test, lfp_list{k}.r_est]));
    axis(ax.lfp(k), 'tight');
    CC = corrcoef(lfp_list{k}.r_test, lfp_list{k}.r_est);
    score.cc.lfp(k) = CC(1,2);
    set(get(get(h(1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off')
    %legend_h = legend('$r(t)$', sprintf('$\\hat{r}(t)$ %.3f',  score.cc.lfp(k)));
    legend_h = legend(sprintf('%.3f',  score.cc.lfp(k)), 'Location', legend_location);
    legend_h.FontSize = legend_FS;
end


%% MUA plots
%set the x labels
set(ax.mua(2:end), 'XTickLabel', '');
%xticklabel = cellfun(@(S) sprintf('%g', 1e3*str2double(S)), ax.mua(1).XTickLabel, 'UniformOutput', false);  % (str) sec --> ms
%ax.mua(1).XTickLabel = xticklabel;  % sec --> ms
xlabel(ax.mua(1), 'Time (ms)');

% set the y labels
set(ax.mua(2:end), 'YTickLabel', '');
%ax.mua(1).YTick = ax.mua(1).YTick([1, end]);  % (Hz)
%ylabel(ax.mua(1), 'Frequency (kHz)');
ylabel(ax.mua(1), 'Z-score');

% title
title(ax.mua(end), 'MUA');
linkaxes(ax.mua);
% axis(ax.mua(end), 'tight');


%% PSTH plots
% get the nonempty\visible PSTH subplots
visible_psth_axes = find( contains(get(ax.psth, 'Visible'), 'on') );
ax0 = visible_psth_axes(2:end);
ax1 = visible_psth_axes(1);


%set the x labels
set(ax.psth( ax0 ), 'XTickLabel', '');
%xticklabel = cellfun(@(S) sprintf('%g', 1e3*str2double(S)), ax.psth(ax1).XTickLabel, 'UniformOutput', false);  % (str) sec --> ms
%ax.psth(ax1).XTickLabel = xticklabel;  % sec --> ms
%xlabel(ax.psth(ax1), 'Lags (ms)');

% set the y labels
set(ax.psth( ax0 ), 'YTickLabel', '');
%ax.psth(ax1).YTick = ax.psth(ax1).YTick([1, end]);  % (Hz)
%ylabel(ax.psth(ax1), 'Frequency (kHz)');

% title
title_psth = title(ax.psth(end), 'Spikes');
% set(findall(ax.psth(1), 'type', 'text'), 'visible', 'on');
set(title_psth, 'visible', 'on');

linkaxes(ax.psth([ax0; ax1]));
% axis(ax.psth(end), 'tight');




%% LFP plots
%set the x labels
set(ax.lfp(2:end), 'XTickLabel', '');
% xticklabel = cellfun(@(S) sprintf('%g', 1e3*str2double(S)), ax.lfp(1).XTickLabel, 'UniformOutput', false);  % (str) sec --> ms
% ax.lfp(1).XTickLabel = xticklabel;  % sec --> ms
%xlabel(ax.lfp(1), 'Lags (ms)');

% set the y labels
set(ax.lfp(2:end), 'YTickLabel', '');
% ax.lfp(1).YTick = ax.lfp(1).YTick([1, end]);  % (Hz)
%ylabel(ax.lfp(1), 'Frequency (kHz)');

% title
title(ax.lfp(end), 'LFP');
linkaxes(ax.lfp);
% axis(ax.lfp(end), 'tight');







