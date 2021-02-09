function [files] = MotSource_interp(cntr_ITD, amp_ITD, ITDmod, SAMmod, SAMdepth, ITDdepth, phase, varargin)
%%% Creates stimuli with both a time-varying ITD and sinusoidally varying
%%% amplitude.  Here, the time-varying ITD is created by upsampling the
%%% original noise signal and then using samples from that at each sampling
%%% point to produce the delay.  The sinusoidal envelope is applied 
%%% afterwards by multiplying the L and R channels by the envelope.  The 
%%% stimulus is played with soundsc before it is saved.
%%% Inputs:
%%%   - cntr_ITD = array of center ITDs for the motion (us)
%%%   - amp_ITD = amplitude of the ITD motion (us).  If it's negative, the
%%%   motion starts at the more negative ITD (closer to ipsi).  If it's 0,
%%%   the noise remains at cntr_ITD
%%%   - mod_rates = array of modulation rates (Hz)
%%%   - SAMmod = modulation frequency for the SAM envelope (Hz)
%%%   - SAMdepth = modulation depth of the SAM envelope (dB, 10log_10)
%%%   - ITDdepth = "modulation depth" of the moving ITD, in dB, using a
%%%   mixing algorithm similar to Siveke et al (2007) and Dietz et al
%%%   (2008)
%%%   - phase = phase difference between the SAM modulation and the DynITD
%%%   modulation (cycles)
%%%
%%% NOTE: The upsampling is done on segments of the original sound in order
%%% to avoid memory errors for long duration sound files.  Once all of the
%%% samples within one segment are retrieved, the next segment is used.

%%% EDIT (7/28/14): Should alternate between ipsi and contra stimuli for varying ITD.
%%% This should be set so the delay is always decreasing for the side with
%%% the varying delay. (+ ITD is toward contra, - ITD = away from contra)
%%% EDIT (7/28/14): IP is in first column, CO is in second.  For
%%% psychophysics, a positive ITD is to the right, negative is to the left.
%%% This way, when the stimuli are played directly to the listener, the
%%% positive and negative ITDs are presented with the proper orientation.
%%% EDIT (11/4/14): Vary delay in the channel where the delay is
%%% decreasing.  Also include low-pass rectangle filters before and after the noise.

%% Initial values
Fs = 100000; % sampling rate (Hz)
tok = randi(100);
both_delay = 0; % flag to specify if both channels have time-varying delay

t_end = 1; % length of time for the stimulus (s)
sect_size = 1000; % number of samples per section over which to interpolate the noise
upsamp = 10; % magnitude by which to upsample
normalize = 1;
flt_noise = 0; % flag to specify if noise should be filtered before interpolation
cutoff_noise = 20000; % cutoff frequency for the low-pass filtered noise, in Hz

files = {};
width = 0.5;

% step_ITD = round(1/Fs*10^6); % ITD step size during ITD motion (us)
motion_type = 'sinusoidal';
AM_type = 'sinusoidal';

if ~isempty(varargin),
    for n = 2:2:length(varargin),
        eval([varargin{n-1} '=varargin{n};']);
    end
end

% Compute the number of sections, based on the section size
% Use 1/2 overlapping sections
% num_sects = 2*t_end*Fs/sect_size-1;
num_sects = 2*t_end*Fs/sect_size;
if floor(num_sects)~=num_sects,
    rem_size = t_end*Fs-floor(num_sects+1)*sect_size/2; % remainder number of samples
    disp(['Remainder samples: ' num2str(rem_size)]);
end

% Check for the smallest ITD increment of the set of modulation rates
t = 0:1/Fs:t_end-1/Fs; % time array
% Determine the shape of the motion
if ITDmod,
    if strcmp(motion_type,'triangle')
        motion = sawtooth(2*pi*ITDmod*t + (0.5+phase)*(2*pi),width);
    elseif strcmp(motion_type,'sinusoidal')
        motion = cos(2*pi*ITDmod*t + phase*(2*pi));
    end
else
    motion = zeros(1,length(t));
end
ITD = -amp_ITD*motion + cntr_ITD;
dITD = diff(ITD);
% Copy the last point for dITD end, in order to make it the same length as
% ITD
dITD = [dITD dITD(end)];
% min_dITD = min(abs(dITD)); % minimum ITD step, in us

% Use the minimum ITD step to determine the upsampling magnitude
% interp_Fs = 10^6/min_dITD;
% upsamp = ceil(interp_Fs/Fs);
disp(['Upsampling magnitude set to: ' num2str(upsamp)]);

    %% Calculate the trajectory of motion
    %motion = sawtooth(2*pi*r*t,0.5);
%         ITD = amp_ITD*motion + cntr_ITD; % convert to ITDs
        
        
%         for tok = 1:1
%% Create the sounds
rng(tok);
bb_noise = randn(1,length(t));
%             CO_chan = bb_noise; % keep the left channel constant
% % Filter the noise at the specified cutoff frequency
if flt_noise,
    BB = fft(bb_noise);
    zero_inds = cutoff_noise*t_end:length(BB)-cutoff_noise*t_end;
    BB(zero_inds) = 0;
    bb_noise = real(ifft(BB));
end

CO_chan = zeros(1,length(bb_noise)); %%% EDITED: 7/26/14
IP_chan = zeros(1,length(bb_noise));

% Iterate through each section

times = zeros(num_sects,1);
for i = 1:floor(num_sects);
    tic
    % Identify the start and end indices of a section
    start_ind = sect_size/2*(i-1)+1;
    end_ind = sect_size+start_ind-1;
    if end_ind > length(bb_noise), %%% EDITED (10/22/14)
        end_ind = length(bb_noise);
    end
    % Create the interpolated noise segment
    % Need to include points before the start_ind of the section, in case
    % the ITD is negative at the start and requires a sampling point before
    % start_ind
    ITDind = ceil(Fs*max(abs(ITD/10^6))); % the number of sample points 
        % equal to the maximal delay imparted by the ITD trajectory
    if (start_ind-ITDind < 0),
        bb_interp = rect_upsample([bb_noise(end-ITDind:end) bb_noise(start_ind:end_ind+ITDind)]',upsamp);
    elseif  (end_ind+ITDind > length(bb_noise)),
%     elseif  (end_ind+ITDind > sect_size),
        bb_interp = rect_upsample([bb_noise(start_ind-ITDind:end_ind) bb_noise(1:ITDind)]',upsamp);
    else
        bb_interp = rect_upsample(bb_noise(start_ind-ITDind:end_ind+ITDind)',upsamp);
    end
    upsind = ITDind*upsamp; % index shift amount, to account for added signal of duration ITDind
    for n = 1:sect_size/2,
        % If the both delay is chosen, and the ITD is increasing (delay
        % decreasing in the contralateral ear), vary the delay in CO_chan
        if both_delay && (dITD(start_ind+n-1)>0),
            IP_chan(start_ind+n-1) = bb_noise(start_ind+n-1);
            CO_chan(start_ind+n-1) = bb_interp(upsind + (n-1)*upsamp+1 + round(ITD(start_ind+n-1)/10^6*Fs*upsamp));
        else
            CO_chan(start_ind+n-1) = bb_noise(start_ind+n-1);
            % Pick the delayed point from the upsampled signal:
            % - upsind (ITDind*upsamp) -- number of samples in the upsampled
            % segment equal to the maximal amount of delay in the ITD
            % trajectory
            % - (n-1)*upsamp+1 -- The sample point in the upsampled signal equal to
            % the sample point in the original noise at n
            % - -ITD(start_ind+n-1,r)/10^6*Fs*upsamp -- the number of samples to
            % move in the upsampled signal, either forward or backward, in
            % order to retrieve the delayed point equivalent to the ITD. Round
            % to the closest sample
            IP_chan(start_ind+n-1) = bb_interp(upsind + (n-1)*upsamp+1 + round(-ITD(start_ind+n-1)/10^6*Fs*upsamp));
        end
    end
    times(i) = toc;
end

disp(['Average elapsed time: ' num2str(mean(times)) ' +/- ' num2str(std(times)) ' s']);
disp(['Total time: ' num2str(sum(times))]);
    
% If the noise should be filtered afterwards
if flt_noise,
    s = [IP_chan' CO_chan'];
    S = fft(s);
    zero_inds = cutoff_noise*t_end:length(S(:,1))-cutoff_noise*t_end;
    S(zero_inds,:) = 0;
    s = real(ifft(S));
    IP_chan = s(:,1)';
    CO_chan = s(:,2)';
end

%% Mix with uncorrelated noise, based on ITDdepth
uncorr = randn(2,length(t));
mix_rat = 1/(1+sqrt(1/(10^(ITDdepth/20))-1)); % compute mixing ratio
CO_chan = CO_chan*mix_rat + uncorr(1,:)*(1-mix_rat);
IP_chan = IP_chan*mix_rat + uncorr(2,:)*(1-mix_rat);

%% Apply amplitude modulation
%         env = (1 + 10^(SAMdepth/10)*cos(2*pi*t*SAMmod + phase*(2*pi)))/2;
if strcmp(AM_type,'sinusoidal'),
%     env = (1 + 10^(SAMdepth/20)*cos(2*pi*t*SAMmod + 0.5*(2*pi)))/2;
    env = 1 + 10^(SAMdepth/20)*cos(2*pi*t*SAMmod + 0.5*(2*pi)); %%% EDITED 9/30/14
        %%% with no SAM, env divides the signal by 1/2
elseif strcmp(AM_type,'triangle'),
    env = (1 + 10^(SAMdepth/20)*sawtooth(2*pi*t*SAMmod,width))/2;
end
CO_chan = CO_chan.*env;
IP_chan = IP_chan.*env;       

%% Save the sounds
if normalize,
    mu = 1/max(max(abs([CO_chan; IP_chan])));

    % Since Impale wants 1 Vpp sounds, need to scale each separately to just
    % under 1 Vpp.  The scaling factors (mu) will then be attached to the wav
    % file and loaded in Impale.  Impale then recognizes that each of the
    % waveforms were scaled differently to reach the same level, and will use
    % mu to scale the waveforms correctly.
    %         CO_EPLchan = CO_chan*mu_CO; % don't need silence, Impale can include it
    %         IP_EPLchan = IP_chan*mu_IP;
    CO_EPLchan = CO_chan*mu;
    IP_EPLchan = IP_chan*mu;
    files = [IP_EPLchan; CO_EPLchan]';
else
    files = [IP_chan; CO_chan]';
end