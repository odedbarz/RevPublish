%
% plot_introduction.m
%
% Description:
% Compares spectrogram reconstructions (estimations) with other DRR conditions.
%
%

clc
fignum = 10;
verbose = 1;

addpath('../../');
setup_environment;

analyze_setup;



%% Plot properties
fontsize = 32;	%32;
fontsize_big = fix(1.5*fontsize);  %42;
fontsize_bigger = fix(2*fontsize); %64;

markersize = 36;
linewidth = 5;


sp = 10;
% idx_sp = sp == splits.idx;      % indices; time indices for speaker SP
idx_sp = idx_fun(sp);      % indices; time indices for speaker SP

if contains(tbl_metadata.fn{sp}, '_F')
    sp_sex = 'Female';
else
   sp_sex = 'Male'; 
end

aux.cprintf('String', '\nInformation about these spectrograms & reconstructions:\n');
fprintf('- Speaker #  : %d\n', sp);
fprintf('- Speaker Sex: %s\n', sp_sex);
fprintf('- Utterance  : %s\n\n', tbl_metadata.txt{sp});

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
% Sest = data.Sest;

clim = [40, 100];   % (dB) color limits

spec_st.f(end) = 8e3;   % change from 7.9931e+03 to 8e3 for graphic clarity 


%%
figh = figure(fignum);
clf;
ax = aux.tight_subplot(2 , 2, [.01 .03], [.12 .08], [.1 .08]);

hz = 1e-3;

% Converting to the sampling rate of the stimulus
n0_16kHz = floor(t0*stim_st.fs) + 1;
n1_16kHz = floor(t1*stim_st.fs) + 1;

tstim = stim_st.t(n0_16kHz:n1_16kHz);
tstim = tstim - tstim(1);
Ystim = @(nn) stim_st.Y(n0_16kHz:n1_16kHz, drr.ordered(nn));
env1  = envelope(Ystim(1), 100, 'rms');
env2  = envelope(Ystim(drr_k), 100, 'rms');
% env   = filter(ones(10,1)/10, 1, abs(ystim));

% [Nfs, Dfs] = rat(4e3/stim_st.fs);
% ystim = resample( stim_st.Y(:,drr.ordered(2)), Nfs, Dfs );


% ax = subplot(2,2,1);
axes(ax(1));
plth = plot(tstim, 0.45*Ystim(1), 'LineWidth', linewidth);
hold on
plot(tstim, env1, 'LineWidth', linewidth);
hold off
set(ax(1), 'FontSize', fontsize);
axis tight
ax(1).XTickLabel = '';
title( 'Dry');
ylabel('Amplitude', 'fontsize', fontsize_big);

% ax(2) = subplot(2,2,2);
axes(ax(2));
plth(2) = plot(tstim, 0.45*Ystim(drr_k), 'LineWidth', linewidth);
hold on
plot(tstim, env2, 'LineWidth', linewidth);
hold off
set(gca, 'FontSize', fontsize);
axis tight
ax(2).XTickLabel = '';
ax(2).YTickLabel = '';
title( sprintf('DRR: %s', drr_k_label) );

% ax(3) = subplot(2,2,3);
axes(ax(3));
[~, surf_h] = spec.plot_spectrogram(t_sp, hz*spec_st.f, Sdry(:,idx_sp),...
    'fontsize', fontsize, 'precision', 2, 'fignum', figh, 'nolabels', 1);
set(gca, 'FontSize', fontsize);
ylabel(ax(3), 'Frequency (kHz)', 'fontsize', fontsize_big);
xlabel(ax(3), 'Time (sec)', 'fontsize', fontsize_big);
caxis(ax(3), clim);    % color boundaries


% ax(4) = subplot(2,2,4);
axes(ax(4));
[~, surf_h(2)] = spec.plot_spectrogram(t_sp, hz*spec_st.f, Sdrr(:,idx_sp),...
    'fontsize', fontsize, 'precision', 2, 'fignum', figh, 'nolabels', 1);
set(gca, 'FontSize', fontsize);
colorbar(ax(4));
ylabel(ax(4), '', 'fontsize', fontsize_big);
ax(4).YTickLabel = '';
xlabel(ax(4), 'Time (sec)', 'fontsize', fontsize_big);
ax(2).YTickLabel = '';
caxis(ax(4), clim);    % color boundaries


linkaxes(ax(1:2));
linkaxes(ax, 'x');
drawnow;
ax(4).Position([1,3]) = ax(2).Position([1,3]);


set(ax, 'FontName', 'Ariel');











