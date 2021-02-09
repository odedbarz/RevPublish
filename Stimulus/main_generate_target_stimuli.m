function trg = main_generate_target_stimuli(Fs, gendre_str, duration_sec, fignum)
%
%   function main_generate_target_stimuli(Fs)
%

if 1 > nargin || isempty(Fs)
    Fs = 100e3;     % (Hz)
end

if 3 > nargin || isempty(duration_sec)
    duration_sec = 3;     % (sec)
end

if 4 > nargin || isempty(fignum)
    fignum = [];     
end

len_y = ceil(Fs * duration_sec);     % (samples)

% Target
rand_opt = 1;   % speech with random phase
% rand_opt = 2;   % speech with shuffled phase

trg.path = './.database/Project1/Spch_(36)sec/';

% Use ONE utterance
switch gendre_str
    case 'F'
        trg.fn = 'TIMIT_train_dr1_FCJF0_SA1.wav'; 
        %trg.fn = '.\.database\TIMIT_Project2\TIMIT_test_dr1_FAKS0_SA2_Fs(100k)Hz.wav'; 
    case 'M'
        trg.fn = 'TIMIT corpus\TIMIT_Males_Fs(100kHz)\TIMIT_train_dr3_MWGR0_SA1.wav';
        %trg.fn = '.\.database\TIMIT_Project2\TIMIT_test_dr1_MSTK0_SA2_Fs(100k)Hz.wav';
end
trg.gendre_str = gendre_str;
[trg.y, Fs_load] = audioread( [trg.path, trg.fn] );

%assert(Fs_load == Fs);
if Fs_load ~= Fs
    [num, den] = rat(Fs/Fs_load);
    trg.y = resample(trg.y, num, den);
end
trg.Fs = Fs;


%% Highpass 
% Perform the HPF on lower (16k Hz) sampling rate
%{
if Fs ~=16e3
    [num, den] = rat(16e3/Fs);
    y16k = resample(trg.y, num, den);
else
    y16k = trg.y;
end

% HPF
trg.hpf.Fs    = 16e3;       % (Hz) 
trg.hpf.Fstop = 50;         % Stopband Frequency
trg.hpf.Fpass = 200;        % Passband Frequency
trg.hpf.Astop = 20;         % Stopband Attenuation (dB)
trg.hpf.Apass = 1;          % Passband Ripple (dB)
[~, trg.hpf.b, trg.hpf.a] = filt.HPF_BW(trg.hpf.Fs, trg.hpf.Fstop, trg.hpf.Fpass,trg.hpf. Astop, trg.hpf.Apass);
y16k_hpf = filtfilt(trg.hpf.b, trg.hpf.a, y16k);

% Back to sampling rate of the signal
trg.y = resample(y16k_hpf, den, num);
%}


%%

% Make sure that the signal <trg.y> is at the desired length
if length(trg.y) > len_y
    fprintf('--> [main_generate_target_stimuli.m]: I had to truncate the signal to the requested length!');
    trg.y = trg.y(1:len_y);
else
    % Padd with zeros
    trg.y = [trg.y; zeros(len_y-length(trg.y),1)];    
end

% Speech with randomized phase
trg.yrphs = shaped_speech_noise( trg.y, rand_opt );

% Calc. the envelope
trg.foct = 1500;     % (Hz) the center frequency of the octave-band filter
[trg.yenv, ~] = filter_speech_env(trg.y, 'foct', trg.foct);


% Modulation (impose an envelope on the speech-noise)
trg.yspch = trg.yenv .* trg.yrphs;




%% Brute force filtering in the frequency domain
trg.hpf.type = 'Brute force';
trg.hpf.Fstop = 30;     % (Hz)
Y = fft(trg.yspch);
W = ones(len_y,1);
win_idx = round(trg.hpf.Fstop/Fs*len_y);
W(1:win_idx) = 0;
W((end-win_idx):end) = 0;
trg.yspch = ifft(Y.*W, 'symmetric');



%% Plot
if isempty(fignum), return; end


tt = linspace(0, 1/Fs*len_y, len_y);
ff = linspace(0, Fs, len_y);

figure(15);
clf;
subplot(2,1,1);
plot(tt, [trg.y, trg.yenv, trg.yspch]);
legend('$y(t)$', sprintf('$y_{env}(t)$: BPF (octave-band, %.1f kHz) + env', 1e-3*trg.foct), '$y_{spch}(t)$: speech-modulated speech-shaped noise');
xlabel('Time (ms)');
ylabel('Amp.');
title('Target (one utterance speech-shaped noise)');
%
subplot(2,1,2);
plot(ff, 20*log10(abs(fft([trg.y, trg.yenv, trg.yspch]))));
xlabel('Frequency (Hz)');
ylabel('Spectrum (dB)');
xlim([0, 16e3/2]);