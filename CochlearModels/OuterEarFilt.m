function Hd = OuterEarFilt(Fs)
%
% Pre-emphasis filter - middle ear filter, 
%   taken from Meddis & O'Mard, 1997
%

% MATLAB Code
% Generated by MATLAB(R) 8.2 and the Signal Processing Toolbox 6.20.
% Generated on: 28-Jul-2014 15:34:35

% Butterworth Bandpass filter designed using FDESIGN.BANDPASS.

% All frequency values are normalized to 1.

N   = 2;       % Order
% Fc1 = 0.0281;  % First Cutoff Frequency
% Fc2 = 0.5313;  % Second Cutoff Frequency
Fc1 = 450;  % First Cutoff Frequency

% Fc2 = 8500;  % Second Cutoff Frequency
Fc2 = 8000;  % Second Cutoff Frequency

% Construct an FDESIGN object and call its BUTTER method.
h  = fdesign.bandpass('N,F3dB1,F3dB2', N, Fc1, Fc2, Fs);
Hd = design(h, 'butter');

% [EOF]
