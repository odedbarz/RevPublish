% main_amgrams.m

clc
verbose = 1;

% Load data into workspace
% analyze_setup;
data = load('D:\Dropbox (Partners HealthCare)\codeOnCloud\RevPublish\_data\Analysis\analyzed_MUA_CC_(04-Nov-2021)_unit_list(103).mat');


%% Prepare
figure(1);
t = data.spec_st.t;
f = data.spec_st.f;

Sstim_dry = data.spec_st.Sft{3};
Sstim_drr = data.spec_st.Sft{4};

dt = diff(data.spec_st.t(1:2));     % temporal sampling rate of the spectrogram

spec.plot_spectrogram(t, f, Sstim_dry);
% title(sprintf('Label: %s\\%%', data.spec_st.labels{k}(1:end-1)));

drr = get_DRR_list_and_indices; 


%% Calculate AMs for all frequency bands
am_nf = 25;
nf = size(Sstim_dry,1);
n_splits = data.splits.n_splits;

Sam_dry = zeros(n_splits, nf, am_nf);
Sam_drr = zeros(n_splits, nf, am_nf);
Xam_dry = zeros(n_splits, nf, am_nf);
Xam_drr = zeros(n_splits, nf, am_nf);

drr_2_use = 5
drr_label   = drr.labels{drr.ordered(drr_2_use)};

Sdry = data.Sest(:,:,1);
Sdrr = data.Sest(:,:,drr_2_use);

for sp = 1:n_splits
    split_idx = sp == data.splits.idx;
    for k = 1:size(Sstim_dry,1)
        [Sam_dry(sp,k,:), ~, am_f] = am_one_freq_band(Sstim_dry(k,~split_idx), dt, am_nf);
        Sam_drr(sp,k,:)            = am_one_freq_band(Sstim_drr(k,~split_idx), dt, am_nf);
        Xam_dry(sp,k,:)            = am_one_freq_band(Sdry(k,~split_idx), dt, am_nf);
        Xam_drr(sp,k,:)            = am_one_freq_band(Sdrr(k,~split_idx), dt, am_nf);    
    end
end

Sam_dry_mean = squeeze(mean(Sam_dry,1));
Sam_dry_std  = squeeze(std (Sam_dry,1));
Sam_drr_mean = squeeze(mean(Sam_drr,1));
Sam_drr_std  = squeeze(std (Sam_drr,1));
Xam_dry_mean = squeeze(mean(Xam_dry,1));
Xam_dry_std  = squeeze(std (Xam_dry,1));
Xam_drr_mean = squeeze(mean(Xam_drr,1));
Xam_drr_std  = squeeze(std (Xam_drr,1));



%% Create the AM-gram for one selected frequency
am_nf = 200;
fband = 16;     % frequency band to analyze/plot

am_stim_avg_dry = squeeze(mean(Sam_dry(:,fband,:),1));   % estimated averaged DRY
am_est_avg_dry = squeeze(mean(Xam_dry(:,fband,:),1));   % estimated averaged DRY
am_est_avg_drr = squeeze(mean(Xam_drr(:,fband,:),1));   % estimated averaged DRR
am_est_std_dry = squeeze(std(Xam_dry(:,fband,:),1));   % estimated averaged DRY
am_est_std_drr = squeeze(std(Xam_drr(:,fband,:),1));   % estimated averaged DRR


figure(10);
errorbar([am_f, am_f], [am_est_avg_dry, am_est_avg_drr], [am_est_std_dry, am_est_std_drr]);
% plot(am_f, [am_est_avg_dry, am_est_avg_drr]);
xlabel('AM Frequency (Hz)');
legend('Reconstructed Dry AMs', sprintf('Reconstructed Reverberant AMs (DRR: %s)', drr_label));
hold on
plot(am_f, am_stim_avg_dry, ':k', 'displayname', 'Avg. Dry-Stimulus AMs')
hold off
title(sprintf('FFT of AMs (frequency band: %.1f Hz)', ...
    data.spec_st.f(fband)));




%%
%{
figure(20);

ax = subplot(1,3,1);
% y1 = max(-100,Xam_dry_mean-Sam_dry_mean);
y1 = Xam_dry_mean./Sam_dry_mean;
spec.plot_spectrogram(am_f, f, y1, 'ax', gca);
xlabel('AM Frequency (Hz)');
colorbar;
title('$\hat{S}_{DRY}/S_{DRY}$');
caxis([0.6, 2.2]);

ax(2) = subplot(1,3,2);
% y2 = max(-100,Xam_drr_mean-Sam_drr_mean);
y2 = Xam_drr_mean./Sam_drr_mean;
spec.plot_spectrogram(am_f, f, y2, 'ax', gca);
xlabel('AM Frequency (Hz)');
colorbar;
title(sprintf('$\\hat{S}_{%s}/S_{%s}$', drr_label, drr_label));
caxis([0.6, 2.2]);

ax(3) = subplot(1,3,3);
% y3 = y1 - y2;
y3 = y2./y1;
spec.plot_spectrogram(am_f, f, y3, 'ax', gca);
xlabel('AM Frequency (Hz)');
colorbar;
title(sprintf('$\\frac{\\hat{S}_{%s}/S_{%s}}{\\hat{S}_{DRY}/S_{DRY}}$',...
    drr_label, drr_label));
caxis([0.6, 2.2]);

linkaxes(ax);
%}



%%
figure(20);
clf;

ax = subplot(1,3,1);
% y1 = max(-100,Xam_dry_mean-Sam_dry_mean);
y1 = Xam_drr_mean./Xam_dry_mean;
spec.plot_spectrogram(am_f, f, y1, 'ax', gca);
xlabel('AM Frequency (Hz)');
colorbar;
title(sprintf('$RMD_{%s}/RMD_{DRY}$', drr_label));

caxis([0.6, 2.2]);

ax(2) = subplot(1,3,2);
% y2 = max(-100,Xam_drr_mean-Sam_drr_mean);
y2 = Sam_drr_mean./Sam_dry_mean;
spec.plot_spectrogram(am_f, f, y2, 'ax', gca);
xlabel('AM Frequency (Hz)');
colorbar;
title(sprintf('$m_{%s}/m_{DRY}$', drr_label));
caxis([0.6, 2.2]);

ax(3) = subplot(1,3,3);
y3 = y1./y2;
spec.plot_spectrogram(am_f, f, y3, 'ax', gca);
xlabel('AM Frequency (Hz)');
colorbar;
title(sprintf('$\\frac {RMD_{%s}/RMD_{DRY}} {m_{%s}/m_{DRY}} $',...
    drr_label, drr_label));
caxis([0.6, 2.2]);

linkaxes(ax);


%%
figure(25);
clf;
freqk = 20;

plot(am_f, [y1(freqk,:)', y2(freqk,:)']);
xlabel('AM Frequency (Hz)');
title(sprintf('AM Ratios (frequency band: %.1f Hz)', ...
    data.spec_st.f(freqk)));
h = legend(sprintf('$RMD_{%s}/RMD_{DRY}$', drr_label),...
    sprintf('$m_{%s}/m_{DRY}$', drr_label));
h.FontSize = 20;
h.Location = 'best';





%%
freqk = 16;
xdry = Xam_dry_mean(freqk,:)' - Sam_dry_mean(freqk,:)';
xdry_std = Xam_dry_std(freqk,:)' - Sam_dry_std(freqk,:)';

xdrr = Xam_drr_mean(freqk,:)' - Sam_drr_mean(freqk,:)';
xdrr_std = Xam_drr_std(freqk,:)' - Sam_drr_std(freqk,:)';

figure(30);
plot(am_f, [xdry, xdrr]);
% errorplot([xdry, xdrr], );

h1 = plot(am_f, [xdry, xdrr]);
hold on
h2 = errorbar([am_f, am_f], [xdry, xdrr], [xdry_std, xdrr_std]);
h2(1).Color = h1(1).Color;
h2(2).Color = h1(2).Color;
% plot(am_f, Sam_dry_mean(freqk,:)', 'k')
hold off



