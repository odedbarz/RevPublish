%
% plot_STRF_array_script.m
% 
% Description:
% Plots array of STRFs for the MUA & spikes
%
% NOTE: lowest subplot is the first channel (Enum==1, i.e., the deepest) 
%


clear ax
ax.mua = [];
ax.lfp = [];
ax.psth = [];
len_mua = length(mua_list);

for k = 1:len_mua
    % MUA
    ax.mua = [ax.mua, subplot(len_mua, 3, 3*len_mua-1 - (3*k-2))];  
    strfpkg.plot_strf(mua_list{k}.strf, spec_st.f, binwidth);
    xlabel('');
    ylabel('');
    
    % LFP
    ax.lfp = [ax.lfp, subplot(len_mua, 3, 3*len_mua+1 - (3*k-2))];
    strfpkg.plot_strf(lfp_list{k}.strf, spec_st.f, binwidth);
    xlabel('');
    ylabel('');
    
    % Spikes (PSTH)
    ax.psth = [ax.psth, subplot(len_mua, 3, 3*len_mua - (3*k-2))];
    if ~isempty(psth_list{k})
        strfpkg.plot_strf(psth_list{k}.strf, spec_st.f, binwidth);
        xlabel('');
        ylabel('');
    else
        set(ax.psth(k), 'Visible', 'Off');
    end
    
end


%% MUA plots
%set the x labels
set(ax.mua(2:end), 'XTickLabel', '');
xticklabel = cellfun(@(S) sprintf('%g', 1e3*str2double(S)), ax.mua(1).XTickLabel, 'UniformOutput', false);  % (str) sec --> ms
ax.mua(1).XTickLabel = xticklabel;  % sec --> ms
xlabel(ax.mua(1), 'Lags (ms)');

% set the y labels
set(ax.mua(2:end), 'YTickLabel', '');
ax.mua(1).YTick = ax.mua(1).YTick([1, end]);  % (Hz)
ylabel(ax.mua(1), 'Frequency (kHz)');

% title
title(ax.mua(end), '$STRF_{MUA}$');
linkaxes(ax.mua);


%% LFP plots
%set the x labels
set(ax.lfp(2:end), 'XTickLabel', '');
xticklabel = cellfun(@(S) sprintf('%g', 1e3*str2double(S)), ax.lfp(1).XTickLabel, 'UniformOutput', false);  % (str) sec --> ms
ax.lfp(1).XTickLabel = xticklabel;  % sec --> ms
%xlabel(ax.lfp(1), 'Lags (ms)');

% set the y labels
set(ax.lfp(2:end), 'YTickLabel', '');
ax.lfp(1).YTick = ax.lfp(1).YTick([1, end]);  % (Hz)
%ylabel(ax.lfp(1), 'Frequency (kHz)');

% title
title(ax.lfp(end), '$STRF_{LFP}$');
linkaxes(ax.lfp);



%% PSTH plots
% get the nonempty\visible PSTH subplots
visible_psth_axes = find( contains(get(ax.psth, 'Visible'), 'on') );
ax0 = visible_psth_axes(2:end);
ax1 = visible_psth_axes(1);


%set the x labels
set(ax.psth( ax0 ), 'XTickLabel', '');
xticklabel = cellfun(@(S) sprintf('%g', 1e3*str2double(S)), ax.psth(ax1).XTickLabel, 'UniformOutput', false);  % (str) sec --> ms
ax.psth(ax1).XTickLabel = xticklabel;  % sec --> ms
%xlabel(ax.psth(ax1), 'Lags (ms)');

% set the y labels
set(ax.psth( ax0 ), 'YTickLabel', '');
ax.psth(ax1).YTick = ax.psth(ax1).YTick([1, end]);  % (Hz)
%ylabel(ax.psth(ax1), 'Frequency (kHz)');

% title
title_psth = title(ax.psth(end), '$STRF_{Spikes}$');
% set(findall(ax.psth(1), 'type', 'text'), 'visible', 'on');
set(title_psth, 'visible', 'on');

linkaxes(ax.psth([ax0; ax1]));












