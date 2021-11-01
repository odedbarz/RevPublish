% main_amgrams.m

clc
verbose = 1;

% Load data into workspace
analyze_setup;


%% Prepare
figure(1);
t = data.spec_st.t;
f = data.spec_st.f;

Sstim_dry = data.spec_st.Sft{3};
Sstim_drr = data.spec_st.Sft{4};

dt = diff(data.spec_st.t(1:2));     % temporal sampling rate of the spectrogram

spec.plot_spectrogram(t, f, Sstim_dry);
% title(sprintf('Label: %s\\%%', data.spec_st.labels{k}(1:end-1)));



%% Create the AM-gram for one selected frequency
am_nf = 200;
fband = 26;
[am_dry, ~, am_f] = am_one_freq_band(Sstim_dry(fband,:), dt, am_nf);

%
drr_2_use = 5
Sdrr = data.Sest(:,:,drr_2_use);
am_est = am_one_freq_band(Sdrr(fband,:), dt, am_nf);

%
figure(10);
semilogy(am_f, [am_dry, am_est]);
xlim([0, 0.5*1/dt]);
xlabel('Frequency (Hz)');
legend('Dry-stimulus AM', 'Est AM');



%% Calculate AMs for all frequency bands
nf = size(Sstim_dry,1);
Sam_dry = zeros(nf, am_nf);
Sam_drr = zeros(nf, am_nf);
Xam_dry = zeros(nf, am_nf);
Xam_est = zeros(nf, am_nf);

drr_2_use = 5
Sdry = data.Sest(:,:,1);
Sdrr = data.Sest(:,:,drr_2_use);

for k = 1:size(Sstim_dry,1)
    [Sam_dry(k,:), ~, am_f] = am_one_freq_band(Sstim_dry(k,:), dt, am_nf);
    Sam_drr(k,:)            = am_one_freq_band(Sstim_drr(k,:), dt, am_nf);
    Xam_dry(k,:)            = am_one_freq_band(Sdry(k,:), dt, am_nf);
    Xam_est(k,:)            = am_one_freq_band(Sdrr(k,:), dt, am_nf);    
end

figure(20);
ax = subplot(1,3,1);
% imagesc(log10(Xam_dry));
spec.plot_spectrogram(am_f, f, log10(Xam_dry), 'ax', gca);
xlabel('AM Frequency (Hz)');
colorbar;
xlim([1, 0.5*1/dt]);

ax(2) = subplot(1,3,2);
% imagesc(log10(Xam_est));
spec.plot_spectrogram(am_f, f, log10(Xam_dry)-log10(Sam_dry), 'ax', gca);
xlabel('AM Frequency (Hz)');
colorbar;
xlim([1, 0.5*1/dt]);

ax(3) = subplot(1,3,3);
% imagesc(max(0, log10(Xam_est)-log10(Xam_dry)));
spec.plot_spectrogram(am_f, f, log10(Xam_est)-log10(Sam_drr), 'ax', gca);
xlabel('AM Frequency (Hz)');
colorbar;
xlim([1, 0.5*1/dt]);
linkaxes(ax)





