% Mesgarani_2014_T_to_DRR.m
%
%

clc
fignum = 10;
verbose = 1;

setup_environment('../');



%% Load data
data_type   = 'MUA';       % {'SU', MUA'}
switch data_type
    case 'SU'
        n_units = 103;  %[10, 25, 50, 103,]; 
        
    case 'MUA'
        n_units = 241;   %[10, 25, 50, 103, 150, 241];    

    otherwise
        error('--> Unrecognized DATA_TYPE!');
        
end


fprintf('Loading the CC arrays...\n');
fn_path = load.path_to_data('Analysis');
fn_name = sprintf('analyzed_cc_%s_units(%d).mat', data_type, n_units);
fn_fullfile = fullfile( fn_path, fn_name );
warning off
data = load(fn_fullfile, 'stim_st', 'spec_st');
warning on




%% Loads stimulus RIR
fn_data = fullfile(load.path_to_data('wav_spch_36sec'), 'spch_36_metadata_new.mat');
dummy   = load(fn_data, 'Trir');
tbl_rir = dummy.Trir;

fs      = tbl_rir.Fs(1);
n_smp   = size(tbl_rir.rir{1}, 1);
t       = (0:n_smp-1)'/fs;

row = 6;

rir = tbl_rir.rir{row}(:,1);
% plot(t, db(rir) );
% aux.hline(-60);

fprintf('\nData\n');
% fprintf('- wtype: %g\n', tbl_rir.wtypes(row));
% fprintf('- dist : %g\n', tbl_rir.dist(row));
% fprintf('- drrL : %g\n', tbl_rir.drrL(row));
disp( tbl_rir(6,:) );



%% Trying to compute DRR directly (no noise)
%n_smp       = tbl_rir.num_taps(1);
duration_sec= 1/tbl_rir.Fs(1) * n_smp;          % (sec)
sigma = 0.300;                            % (sec)
% t           = linspace(0, duration_sec, n_smp);
delta = double(1 == 1:length(t))';

Ed = db( sum( delta.^2 ) );     % Direct energy (equals 0)

% exp_fun = @(T,SIG) 1/sqrt(2*pi*SIG^2)*exp(-0.5*(T./SIG).^2);
exp_fun = @(T,SIG) exp(-T./SIG);
rir_exp = exp_fun(t,sigma);
%     rir_exp_noise = rir_exp;    '$$$ DEBUG $$$'
noise = randn(length(t),1);
rir_exp_noise = max(0, noise .* rir_exp);
    
%         '### DEBUG ###'
%         %rir_exp_noise = delta + rir_exp_noise;  '### DEBUG ###'
%         %rir_exp_noise = std(rir_exp_noise)*rir_exp_noise;
%         rir_exp_noise(1) = 0;
%         rir_exp_noise = rir_exp_noise/max(rir_exp_noise);
%         rir_exp_noise = delta + rir_exp_noise;  
%         
% A       = trapz(t, rir_exp_noise);
% Er      = db( A );      % Reverberation energy
% drr_exp = Ed - Er;      % (dB)

rir_exp_noise(1) = 0;
rir_exp_noise = rir_exp_noise/max(rir_exp_noise);
rir_exp_noise = delta + rir_exp_noise;  
drr_exp = Direct_to_Reverberation(rir_exp_noise, 1);

fprintf('\nTrying to compute DRR directly (no noise)\n');
fprintf('- DRR exp: %.2f dB\n', drr_exp);





% Plot RT60
figure(fignum);
clf;
eps_ = 1e-5;

num_taps = tbl_rir.num_taps(1,1);
Fs = tbl_rir.Fs(1,1);

t = linspace(0, num_taps/Fs, num_taps)';
warning off
plth = plot(t, db( eps_+rir ), 'DisplayName', sprintf('RIR, (DRR: %.2f dB; room-image)', tbl_rir.drrL(row)));
warning on
hold on
plth(2) = plot(t, db(rir_exp_noise), 'DisplayName', 'exp RIR $\times$ noise');
plth(3) = plot(t, db(rir_exp), 'DisplayName', 'exp RIR');
plth(4) = plot(t, envelope(db( eps_+rir ), 250, 'peak'), 'DisplayName', 'Envelope of the RIR (room-image)');
hold off
ylabel('RIR (dB)');
xlabel('Time (sec)');

aux.ctitle('Room Impulse Response',...
    sprintf('DRR: %.2f dB,  Gauss-DRR: %.2f dB', tbl_rir.drrL(row), drr_exp));

% xlim([0, t(end)]);
% xlim([0.3, 1.5]);
ylim([-117, -2.4]);

hold on
line_h = line(xlim, -60*[1 1]);
line_h.Color = 'k';
line_h.LineStyle = '--';
line_h.LineWidth = 2;
text_h = text(0.9*max(xlim), -60, '$T_{60}$');
text_h.Interpreter = 'latex'; 
text_h.FontSize = 32; 
text_h.VerticalAlignment = 'bottom'; 
hold off


legend(plth, 'Location', 'southwest');


% Plot the convolution kernels
figure(2+fignum);
clf;
plth = plot(t, rir_exp_noise, 'DisplayName', 'exp RIR $\times$ noise');
warning on
hold on
plth(2) = plot(t, rir_exp, 'DisplayName', 'exp RIR');
plth(3) = plot(t, rir, 'DisplayName', sprintf('RIR, (DRR: %.2f dB; room-image)', tbl_rir.drrL(row)));
hold off
legend(plth, 'Location', 'northeast');
title('Convolution Kernels');
ylabel('Amp.');
xlabel('Time (sec)');
xlim([0.0, 1.0]);






%% Plot the STIMULI + REVERBERATIONs wave forms
y_dry = data.stim_st.Y(:,3);
y_dry = y_dry/max(abs(y_dry));
y_dry_pad = [y_dry; zeros(n_smp,1)];

y_rir = fftfilt( rir, y_dry_pad );
y_rir = y_rir(1:length(y_dry));
y_rir = y_rir/max(abs(y_rir));
y_rir = circshift(y_rir, -34);

y_rir_exp = fftfilt( rir_exp_noise, y_dry_pad );
y_rir_exp = y_rir_exp(1:length(y_dry));
y_rir_exp = y_rir_exp/max(abs(y_rir_exp));

% Plot RT60
figure(10+fignum);
clf;
plth = plot(data.stim_st.t, y_dry, 'DisplayName', 'Dry Stimulus');
warning on
hold on
plth(2) = plot(data.stim_st.t, y_rir_exp, 'DisplayName', 'Stimulus * Gaussian kernel (+ noise)');
plth(3) = plot(data.stim_st.t, y_rir, 'DisplayName', sprintf('Stimulus * room-image-RIR (DRR: %.2f dB)', tbl_rir.drrL(row)));
hold off
h = legend(plth, 'Location', 'northeast');
h.FontSize = 24;
xlim([0.0, 3.0]);
ylabel('Amp.');
xlabel('Time (sec)');




%% Load stimuli & measurement data
fprintf('--> Load stimuli & spectrograms\n');
spectrogram_type = 'gammatone';      % {['matlab'], 'stft', 'multitaper', 'gammatone'}
f_scale     = 'log';	% {['lin'], 'log', 'erb'}
n_bands     = 30;       % (1x1) # of bins along the frequency domain of the spectrogram
binwidth    = 5;        % (ms) binwidth of the resulted spectrogram 
win_size_ms = nan;      % (ms) temporal window size over which to calc the spectrogram; 
                        %      'gammatone' filterbanks do not use it!
lowfreq     = 200;      % (Hz)
highfreq    = 8840;     % (Hz) %8800;     
nw          = [];       % applies only for SPECTROGRAM_TYPE = 'multitaper'
neurons     = 1;        % neuron #1 (#115) is for stimulus of 36 sec (40 sec)
spectral_diff= 0;       % (logical) perform derivative (DIFF) along the frequency domain
hpf_pole    = nan;

% Spectrograms
[Sdry, spec_st] = spec.spectrogram(y_dry, fs, ...
    'n_bands', n_bands,...
    'lowfreq', lowfreq,...
    'highfreq', highfreq,...
    'overlap_ratio', 0.80,...
    'binwidth', binwidth,...
    'win_size_ms', win_size_ms, ...
    'nw', nw,...                only for spectrogram_type== MULTITAPER
    'f_scale', f_scale,...
    'db_floor', -100, ...  % (dB)
    'duration_ms', 1e3*duration_sec,...
    'method', spectrogram_type, ...
    'fignum', [] ...
);

[Srir, spec_st] = spec.spectrogram(y_rir, fs, ...
    'n_bands', n_bands,...
    'lowfreq', lowfreq,...
    'highfreq', highfreq,...
    'overlap_ratio', 0.80,...
    'binwidth', binwidth,...
    'win_size_ms', win_size_ms, ...
    'nw', nw,...                only for spectrogram_type== MULTITAPER
    'f_scale', f_scale,...
    'db_floor', -100, ...  % (dB)
    'duration_ms', 1e3*duration_sec,...
    'method', spectrogram_type, ...
    'fignum', [] ...
);

[Sexp, spec_st] = spec.spectrogram(y_rir_exp, fs, ...
    'n_bands', n_bands,...
    'lowfreq', lowfreq,...
    'highfreq', highfreq,...
    'overlap_ratio', 0.80,...
    'binwidth', binwidth,...
    'win_size_ms', win_size_ms, ...
    'nw', nw,...                only for spectrogram_type== MULTITAPER
    'f_scale', f_scale,...
    'db_floor', -100, ...  % (dB)
    'duration_ms', 1e3*duration_sec,...
    'method', spectrogram_type, ...
    'fignum', [] ...
);


%% Plot DRY spectrogram
figure(20+fignum);
clf;
ax = subplot(3,1,1);
plth = spec.plot_spectrogram(spec_st.t, spec_st.f, Sdry, 20+fignum);
plth(1).XTickLabel = '';
xlabel('');
ylabel('Dry');

ax(2) = subplot(3,1,2);
plth(2) = spec.plot_spectrogram(spec_st.t, spec_st.f, Srir, 20+fignum);
plth(2).XTickLabel = '';
xlabel('');
ylabel( aux.ctitle('Frequency (kHz)', sprintf('DRR: %.2f dB', tbl_rir.drrL(row)) ) );

ax(3) = subplot(3,1,3);
plth(3) = spec.plot_spectrogram(spec_st.t, spec_st.f, Sexp, 20+fignum);
ylabel(sprintf('Gauss-RIR (T: %g ms)', 1e3*sigma));

linkaxes(ax, 'x');
xlim([0.0, 3.0]);
set(ax, 'FontSize', 24);













