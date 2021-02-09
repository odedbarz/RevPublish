% plot_spectrograms.m
%
% Plots stimuli spectrograms and reconstructed spactrograms.


clc
fignum = 10;
verbose = 1;

setup_environment('../');

analyze_setup;



%% Plot properties
fontsize = 32;	%32;
fontsize_big = fix(1.5*fontsize);  %42;
fontsize_bigger = fix(3*fontsize); %64;

markersize = 36;
linewidth = 5;


sp = 10;
% idx_sp = sp == splits.idx;      % indices; time indices for speaker SP
idx_sp = idx_fun(sp);      % indices; time indices for speaker SP

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

ax = zeros(5,2);
surf_h = zeros(5,2);


for rv = 1:n_drr
    rvi = drr.sortby(rv);
    
    Sdrr = squeeze( spec_st.Sft{rvi}(:,idx_sp) );
    Sest = squeeze( data.Sest(:,idx_sp,rv) );
        
    ax(rv,1) = subplot(5,2, 5*2 - (1+2*(rv-1)));    
    [~, surf_h(rv,1)] = spec.plot_spectrogram(t_sp, hz*spec_st.f, Sdrr,...
        'fontsize', fontsize, 'precision', 2, 'fignum', figh, 'nolabels', 1);
    
    ax(rv,2) = subplot(5,2, 5*2-2*(rv-1));
    [~, surf_h(rv,2)] = spec.plot_spectrogram(t_sp, hz*spec_st.f, Sest,...
        'fontsize', fontsize, 'precision', 2, 'fignum', figh, 'nolabels', 1);

    ylabel(ax(rv,1), drr.labels{rvi}, 'FontSize', fontsize_big);
    
end
           

ax = flipud(ax);
arrayfun(@(X) caxis(X, clim), ax);   % (dB) set the range of colors

set(ax, 'FontSize', fontsize);
set(ax(1:4,:), 'XTickLabel', '');
set(ax(:,2), 'YTickLabel', '');
xlabel(ax(end,1), 'Time (sec)', 'FontSize', fontsize_big);
xlabel(ax(end,2), 'Time (sec)', 'FontSize', fontsize_big);

% Add the titles
title(ax(1,1), 'Stimulus', 'FontSize', fontsize_big);
title(ax(1,2), 'Reconstruction', 'FontSize', fontsize_big );

aux.abc(ax(1,:), 'fontsize', fontsize_bigger, 'location', 'northwestoutside');
linkaxes(ax);       
        
set(ax(:,1), 'YTick', log2( hz*[spec_st.f(1), spec_st.f(end)] ) );
yticklabels = arrayfun(@(X) num2str(X,.0), hz*[spec_st.f(1), spec_st.f(end)], 'UniformOutput', false);
set(ax(:,1), 'YTickLabel', yticklabels);

% Add a colorbar and make sure that the axes doesn't change in size (push 
% the colorbar "outside").
if 1
    ax_j = ax(end,end);
    pos = get(ax_j, 'Position');
    
    colorbar(ax_j);
    set(ax_j, 'Position', pos);
end
     





%%   
figh2 = figure(2+fignum);
clf;

ax2 = [];
ylabelh = [];
colors = colormap('jet');

for k = 1:n_drr
    ax2(k) = subplot(n_drr,1,k);
    set(ax2(k), 'XTickLabel', '');
    set(ax2(k), 'YTickLabel', '');
    set(ax2(k), 'Visible', 'Off');
    %set(ax2(k), 'Color', colors(fix(length(colors)/n_drr), :) );
    
    xlim([0, 1]);
    ylim([0, 1]);
    
    ylabelh = ylabel(ax2(k), data.tbl.CC.Row{n_drr-k+1});
    ylabelh.FontSize = round(0.8*fontsize_big);
    ylabelh.Visible = 'On';
    
    
    txt = sprintf('%.2f', data.tbl.CC{n_drr-k+1,sp});
    txth = text(0.5, 0.5, txt);
    txth.FontSize = fix(1.5*fontsize_big);
    txth.HorizontalAlignment = 'center';
    %txth.FontWeight = 'bold';
    %txth.Color = colors(fix(length(colors)/n_drr*k), :);
    
end

% titleh = title(ax2(1), '$CC(S_{dry}$ vs. $\hat{S}_{drr})$', 'FontSize', round(1.2*fontsize_big));
% titleh.Visible = 'On';
figh2.Position(2:4) = [1, 250, figh.Position(4)];










