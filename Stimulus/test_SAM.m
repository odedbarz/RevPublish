%
% test_SAM.m
%
% Create SAM noise signals
%

clc

Fs = 100;               % (kHz)
t_stim = 400;           % (msec) stimulus time
Nt = round(t_stim*Fs);

tt = linspace(0, t_stim, Nt)';
ff = linspace(0, Fs, Nt)';
st = @(f) 1/Nt*sin(2*pi*f*tt);       % [f] = kHz

fm = 0.2;   % (kHz) frequency modulation
stfm = st(fm);

%% Create noise tokens of SAM noise
num_tokens = 50;

tokens = cell2mat(arrayfun(@(X) randn(Nt,1), 1:num_tokens, 'Uniform', 0 ));

sam_mtx = tokens .* stfm;
sam_avg = mean(sam_mtx, 2);


%% Plot
figure(11);
clf;
subplot(1,2,1);
plot(tt, stfm, tt, sam_avg);
xlabel('Time (ms)');
ylabel('s(t)');
legend('s(t)', sprintf('<SAM>_{%d}', num_tokens));
xlim([0, 5/fm]);
subplot(1,2,2);
plot(ff, log10(abs(fft(stfm))), 'DisplayName', '|DFT(s(t))|');
xlim([0, 1.25*fm]);
xlabel('Frequency (kHz)');
ylabel('Mag. log_{10}()');

% figure(11);
subplot(1,2,2);
hold on
plot(ff, log10(abs(fft(sam_avg))), 'DisplayName', sprintf('<SAM>_{%d}', num_tokens) );
hold off
legend('show', 'Location', 'southeast');






















