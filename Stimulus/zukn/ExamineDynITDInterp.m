function ExamineDynITDInterp(s,varargin)
% Plot a spectrogram of the signal s (DynITD interpolated signal)
% In the DynITD stimulus, the ipsilateral side (first channel) has the
% time-varying delay

% Variables
Fs = 100000;
win_size = 30;
av_size = 3;
rms_size = 200;

% Parse varargin
if ~isempty(varargin),
    for n = 2:2:length(varargin),
        eval([varargin{n-1} '=varargin{n};']);
    end
end

sz = size(s);
if min(sz)==1,
    [~,F,T,P] = spectrogram(s,win_size,0,win_size,Fs);
else
    [~,F,T,P] = spectrogram(s(:,1),win_size,0,win_size,Fs);
end
flS = filter2(ones(av_size),P);

% Plot spectrum
figure
imagesc(T,flipud(F),flipud(flS));
xlabel('Time (s)');
ylabel('Frequency (Hz)');

% Plot running rms
avs = filter(ones(rms_size,1),1,s(:,1).^2);
rms = sqrt(avs);
figure
plot((0:size(s,1)-1)/Fs, rms);
xlabel('Time (s)');
ylabel('RMS');