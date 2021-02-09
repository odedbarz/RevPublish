function itr = main_generate_interferer_stimuli(path_to_wav_files, Fs, duration_sec, fignum)
%
%   function main_generate_interferer_stimuli(Fs)
%

if 2 > nargin || isempty(Fs)
    Fs = 100e3;     % (Hz)
end

if 3 > nargin || isempty(duration_sec)
    duration_sec = 3;     % (sec)
end

if 4 > nargin || isempty(fignum)
    fignum = [];     
end

len_y = ceil(Fs * duration_sec);     % (samples)

% files.path   = {'D:\DataBases\TIMIT\Female', 'D:\DataBases\TIMIT\Males'}; 
files.path = path_to_wav_files;
files.all = extract_wav_files( files.path );
assert(~isempty(files.all), '--> ERROR at [main_generate_interferer_stimuli.m] could not find any files!');

[Ywav, fs_wav] = load_wav_files( files.all );

[yavg_fs, num_wav, ~] = averaged_speech_noise( Ywav );
fprintf('--> # of wav files used: %d\n', num_wav);

[num, den] = rat(Fs/fs_wav);
itr.yavg_long = resample(yavg_fs, num, den);
itr.Fs     = Fs;

fprintf('--> Loaded sample rate: %g Hz\n', fs_wav);
fprintf('--> New sample rate   : %g Hz (resample)\n', Fs);


% Make sure that the signal <itr.y> is at the desired length
if length(itr.yavg_long) > len_y
    warning('--> [main_generate_interferer_stimuli.m]: I had to TRUNCATE the signal to the requested length!');
    itr.yavg = itr.yavg_long(1:len_y);
else
    % Padd with zeros
    warning('--> [main_generate_interferer_stimuli.m]: I had to PAD the signal with zeros to the requested length!');    
    itr.yavg = [itr.yavg_long; zeros(len_y-length(itr.yavg_long),1)];    
end






%% Plot
if isempty(fignum), return; end

figure(11);
% clf;
subplot(2,1,1);
plot(tt, itr.yavg);
xlabel('Time (ms)');
ylabel('Amp.');
title('Interferer (one utterance speech-shaped noise)');
%
subplot(2,1,2);
plot(ff, 20*log10(abs(fft(itr.yavg))));
xlabel('Frequency (Hz)');
ylabel('Spectrum (dB)');
xlim([0, 16e3/2]);

