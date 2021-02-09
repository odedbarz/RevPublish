function [files] = MotSource_ITDSweep(cntr_ITD, amp_ITD, mod_rate, varargin)
%%% Creates stimuli with consisting of a single sweep across the ITD range
%%% 
%%% Inputs:
%%%   - cntr_ITD = array of center ITDs for the motion (us)
%%%   - amp_ITD = amplitude of the ITD motion (us).  If it's positive, the
%%%   motion starts at the more negative ITD (closer to ipsi).  If it's 0,
%%%   the noise remains at cntr_ITD
%%%   - mod_rate = array of modulation rates (Hz)

%%% EDITED (8/16/14): The stimulus duration 'dur' is separate from the
%%% sweep duration 't_end'.  If dur~=t_end, then the sweep will be followed
%%% by silence for the rest of the duration.

%% Initial values
Fs = 100000; % sampling rate (Hz)
c = 340.29; % speed of sound in m/s (taken from google search)
% tok = 42;
tok = randi(100);

width = 0.5;
step_ITD = round(1/Fs*10^6); % ITD step size during ITD motion (us)
sweep_type = 'constant';

t_end = 1/(2*mod_rate); % compute duration of sweep, based on the modulation rate specified
dur = 0.25; % duration of the entire sweep, in s (if t_end ~= dur, then the sweep is
    % followed by silence

if ~isempty(varargin),
    for n = 2:2:length(varargin),
        eval([varargin{n-1} '=varargin{n};']);
    end
end

%% Calculate the trajectory of motion
t = 0:1/Fs:t_end-1/Fs;
if mod_rate,
    if strcmp(sweep_type,'constant')
        motion = sawtooth(2*pi*mod_rate*t + (0.5)*(2*pi),width);
    elseif strcmp(sweep_type,'step')
        motion = [ones(1,length(t)/2) -ones(1,length(t)/2)];
    else
        error('Unknown sweep type');
    end
else
    motion = zeros(1,length(t));
end

%         ITD = amp_ITD*motion + cntr_ITD; % convert to ITDs
ITD = -amp_ITD*motion + cntr_ITD;

%% Find the amount of time per bin
%%% Identify the time points in the ITD trajectory between which an
%%% ITD should be set in a stimulus
%%% -- Round each point in the ITD trajectory to the nearest multiple of the step
%%% size.
%%% -- Find the indeces where there is a change in ITD using the
%%% difference of the 'stepped' ITD trajectory.
%%% -- Append 1 to the front of the array, and the last time index
%%% to the end of the array.
%%% To compute the length of the noise between ITD steps:
%%% -- Set the start index to the prior ITD step index + 1, and
%%% the end index to the next ITD step index
ITD_st = round(ITD/step_ITD)*step_ITD;
dITD = diff(ITD_st);
dITD_ind = find(dITD~=0);
dITD_ind = [1 dITD_ind length(t)];

%         for tok = 1:1
%% Create the sounds
rng(tok);
bb_noise = randn(1,length(t));
%             CO_chan = bb_noise; % keep the left channel constant
% CO_chan = zeros(1,length(bb_noise)); %%% EDITED: 7/26/14
% IP_chan = zeros(1,length(bb_noise));
CO_chan = zeros(1,dur*Fs);
IP_chan = zeros(1,dur*Fs);

% Iterate through ITD steps, designated by dITD_ind
for i = 1:length(dITD_ind)-1;
    start_ind = dITD_ind(i)+1;
    end_ind = dITD_ind(i+1);
    ITD_i = ITD_st(start_ind); % get the ITD at the index
        % corresponding to t_incr(i)
    ITD_sf = round(ITD_i/10^6*Fs); % compute the number of indeces for the shift
%                 ITD_sf = ITD_i/10^6*Fs; 
    if (start_ind-ITD_sf < 1) || (end_ind-ITD_sf > length(bb_noise)),
        bb_wrapped = [bb_noise(length(bb_noise)-abs(ITD_sf)+1:end) ...
            bb_noise bb_noise(1:abs(ITD_sf))];
        st_i = start_ind+abs(ITD_sf); % correct for wrapping
        end_i = end_ind+abs(ITD_sf);
    else
        bb_wrapped = bb_noise;
        st_i = start_ind;
        end_i = end_ind;
    end
    if dITD_ind(i) > 0, % contra experiencing decreasing delay
        IP_chan(start_ind:end_ind) = bb_wrapped(st_i-ITD_sf:end_i-ITD_sf);
        CO_chan(start_ind:end_ind) = bb_wrapped(st_i:end_i);
    else % ipsi experiencing decreasing delay
        CO_chan(start_ind:end_ind) = bb_wrapped(st_i-(-ITD_sf):end_i-(-ITD_sf));
        IP_chan(start_ind:end_ind) = bb_wrapped(st_i:end_i);
    end
%                 IP_chan(start_ind:end_ind) = bb_wrapped(st_i-ITD_sf:end_i-ITD_sf);
end      

%% Save the sounds
mu = 1/max(max(abs([CO_chan; IP_chan])));

% Since Impale wants 1 Vpp sounds, need to scale each separately to just
% under 1 Vpp.  The scaling factors (mu) will then be attached to the wav
% file and loaded in Impale.  Impale then recognizes that each of the
% waveforms were scaled differently to reach the same level, and will use
% mu to scale the waveforms correctly.
CO_EPLchan = CO_chan*mu;
IP_EPLchan = IP_chan*mu;

files = [IP_EPLchan; CO_EPLchan]';
