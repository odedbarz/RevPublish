% plot_spectrograms.m
%
% Plots stimuli spectrograms, theri reconstructed spactrograms, and a column of
% text with the respective CCs.


clc
fignum = 10;
verbose = 1;

addpath('../../');
setup_environment;

analyze_setup;



%% Plot properties
fontsize = 32;	%32;
fontsize_big = fix(1.5*fontsize);  %42;
fontsize_bigger = fix(3*fontsize); %64;

markersize = 36;
linewidth = 5;


sp = 10;                % choose a speaker to show
idx_sp = idx_fun(sp); 	% indices; time indices for speaker SP

tidx = t(idx_sp);       % (sec)
t0 = tidx(1);           % (sec)
t1 = tidx(end);         % (sec)
tidx = tidx-t0;         % (sec)

% 1:{'Dry'}, 2:{'9.4 dB'}, 3:{'4.8 dB'}, 4:{'-2.5 dB'}, 5:{'-8.2 dB'}
drr_k  = 5;      
drr_dry_idx = drr.ordered(1);
drr_k_idx = drr.ordered(n_drr);
drr_k_label = drr.labels{drr_k_idx};


fs = 1/(1e-3*spec_st.binwidth);     % (Hz) sampling rate of the spectrogram
t_sp = nonzeros(idx_sp*1/fs);        % (sec) time for the SP speaker
t_sp = t_sp - t_sp(1);

Sdry = spec_st.Sft{drr_dry_idx};
Sdrr = spec_st.Sft{drr_k_idx};
Sest = data.Sest;

% spec_thr    = 20;  % (dB) threshold\remove from spectrogram before plotting it
% Sdry_ = max(0, Sdry - spec_thr);
% Sdrr_ = max(0, Sdrr - spec_thr);

clim = [40, 100];   % (dB) color limits

if contains(tbl_metadata.fn{sp}, '_F')
    sp_sex = 'Female';
else
   sp_sex = 'Male'; 
end



%%
aux.cprintf('String', '\nInformation about these spectrograms & reconstructions:\n');
fprintf('- Speaker #  : %d\n', sp);
fprintf('- Speaker Sex: %s\n', sp_sex);
fprintf('- Utterance  : %s\n\n', tbl_metadata.txt{sp});





%% Plot: spectrograms of stimuli vs. reconstructions
figh = figure(fignum);
clf;


hz = 1e-3;  % (Hz)

ax = zeros(5,3);
surf_h = zeros(5,2);

n_drr = drr.n_drr;
subx = 5;
suby = 24;

for rv = 1:n_drr
    rvi = drr.ordered(end-rv+1);
    
    Sdrr = squeeze( spec_st.Sft{rvi}(:,idx_sp) );
    Sest = squeeze( data.Sest(:,idx_sp,n_drr-rv+1) );       % Sest is already sorted
        
    ax(rv,1) = subplot( subx, suby, suby*(rv-1)+(1:10) );    
    [~, surf_h(rv,1)] = spec.plot_spectrogram(t_sp, hz*spec_st.f, Sdrr,...
        'fontsize', fontsize, 'precision', 2, 'fignum', figh, 'nolabels', 1);
    
    ax(rv,2) = subplot( subx, suby, suby*(rv-1)+(13:22) );
    [~, surf_h(rv,2)] = spec.plot_spectrogram(t_sp, hz*spec_st.f, Sest,...
        'fontsize', fontsize, 'precision', 2, 'fignum', figh, 'nolabels', 1);

    ylabel(ax(rv,1), drr.labels{rvi}, 'FontSize', fontsize_big);
    
    % *** Add the text labels of the CCs ***
    ax(rv,3) = subplot( subx, suby, suby*(rv-1)+(23:24) );
    %set(ax(rv,3), 'XTickLabel', '');
    %set(ax(rv,3), 'YTickLabel', '');
    %set(ax(rv,3), 'Visible', 'Off');
    set(ax(rv,3), 'XColor', 'None');
    set(ax(rv,3), 'YColor', 'None');
    set(ax(rv,3), 'Color', 'None');
    xlim([0, 1]);
    ylim([0, 1]);
    txt = sprintf('%.2f', data.tbl.CC{n_drr-rv+1,sp});
    txth = text(1.0, 0.5, txt);
    txth.FontSize = fix(1.5*fontsize_big);
    txth.HorizontalAlignment = 'center';
    
end
    
arrayfun(@(X) caxis(X, clim), ax(:,1:2));   % (dB) set the range of colors

set(ax(:,1:2), 'FontSize', fontsize);
set(ax(1:4,1:2), 'XTickLabel', '');
set(ax(:,2), 'YTickLabel', '');
xlabel(ax(end,1), 'Time (sec)', 'FontSize', fontsize_big);
xlabel(ax(end,2), 'Time (sec)', 'FontSize', fontsize_big);

% Add the titles
title(ax(1,1), 'Stimulus', 'FontSize', fontsize_big);
title(ax(1,2), 'Reconstruction', 'FontSize', fontsize_big );
% *** for CC TEXT ***
title_h = title(ax(1,3), 'CC', 'FontSize', fontsize_big );    
title_h.Position = [txth.Position(1), title_h.Position(2), title_h.Position(3)];

aux.abc(ax(1,1:2), 'fontsize', fontsize_bigger, 'location', 'northwestoutside');
linkaxes(ax(:,1:2));       
        
set(ax(:,1), 'YTick', log2( hz*[spec_st.f(1), spec_st.f(end)] ) );
yticklabels = arrayfun(@(X) num2str(X,.0), hz*[spec_st.f(1), spec_st.f(end)], 'UniformOutput', false);
set(ax(:,1), 'YTickLabel', yticklabels);

% Add a colorbar and make sure that the axes doesn't change in size (push 
% the colorbar "outside").
if 1
    ax_j = ax(5,1);
    pos = get(ax_j, 'Position');
    
    colorbar(ax_j);
    set(ax_j, 'Position', pos);
end
     









