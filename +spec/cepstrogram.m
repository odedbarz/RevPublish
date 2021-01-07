%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Cepstrogram with MATLAB Implementation     %
%                                                %
% Author: Ph.D. Eng. Hristo Zhivomirov  08/25/16 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [C, q, t] = cepstrogram(x, win, binwidth, fs, min_ceps_ms)
%
% function [C, q, t] = cepstrogram(x, win, binwidth, fs, min_ceps_ms)
%
% OLD: function: [C, q, t] = cepstrogram(x, win, hop, fs)

% Input:
% x - signal in the time domain
% binwidth - (ms) final bin width along the time axis (x-axis)
% win - analysis window function
% hop - hop size
% fs - sampling frequency, Hz
% min_ceps_sec - minimum time for cepstrum coefficient
%
% Output:
% C - real cepstrum-matrix (only unique points, 
%     time across columns, quefrency across rows)
% q - quefrency vector, s
% t - time vector, s

if 5 > nargin
    min_ceps_ms = []; 
end

hop = floor(fs * binwidth*1e-3);  

% representation of the signal as column-vector
x = x(:);

% determination of the signal length 
xlen = length(x);

% determination of the window length
wlen = length(win);

% stft matrix size estimation and preallocation
NUP = ceil((1+wlen)/2);     % calculate the number of unique fft points
L = 1+fix((xlen-wlen)/hop); % calculate the number of signal frames
C = zeros(NUP, L);          % preallocate the stft matrix

% STFT (via time-localized FFT)
for n = 0:L-1
    % windowing
    xw = x(1+n*hop : wlen+n*hop).*win;
    
    % cepstrum calculation
    c = real(ifft(log(abs(fft(xw)))));
    
    % update of the cepstrum-matrix
    C(:, 1+n) = c(1:NUP);
end

% calculation of the time and quefrency vectors
t = (wlen/2:hop:wlen/2+(L-1)*hop)/fs;
q = (0:NUP-1)/fs;   % quefrency


% some conditioning
if ~isempty(min_ceps_ms)
    min_ceps_sec = 1e-3*min_ceps_ms;
    C = C(q > min_ceps_sec, :);   	% ignore all cepstrum coefficients for 
                                  	% quefrencies bellow 0.5 ms  
    q = q(q > min_ceps_sec);      	% ignore all quefrencies bellow 0.5 ms
end

q = q*1000;                         % convert the quefrency to ms


end