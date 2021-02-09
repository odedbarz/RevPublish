function [files] = MotSource_upsampnoise(cntr_ITD, amp_ITD, ITDmod, SAMmod, SAMdepth, ITDdepth, phase, varargin)
%%% Creates stimuli with both a time-varying ITD and sinusoidally varying
%%% amplitude.  Here, the time-varying ITD is created by starting with a
%%% upsampled noise signal and then using samples from that at each sampling
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

%% Initial values
Fs = 100000; % sampling rate (Hz)
tok = randi(100);

t_end = 1; % length of time for the stimulus (s)
sect_size = 1000; % number of samples per section over which to interpolate the noise
upsamp = 10; % magnitude by which to upsample

files = {};
width = 0.5;

% step_ITD = round(1/Fs*10^6); % ITD step size during ITD motion (us)
motion_type = 'triangle';
AM_type = 'sinusoidal';

if ~isempty(varargin),
    for n = 2:2:length(varargin),
        eval([varargin{n-1} '=varargin{n};']);
    end
end

% Compute the number of sections, based on the section size
% Use 1/2 overlapping sections
num_sects = 2*t_end*Fs/sect_size-1;
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
% dITD = diff(ITD);
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
ITDind = ceil(Fs*max(abs(ITD/10^6))); % the number of sample points 
    % equal to the maximal delay imparted by the ITD trajectory
rng(tok);
bb_ups = randn(1,(length(t)+2*ITDind)*upsamp); % may run out of memory with larger upsampling rates
%             CO_chan = bb_noise; % keep the left channel constant
CO_chan = zeros(1,length(t)); %%% EDITED: 7/26/14
IP_chan = zeros(1,length(t));

% Iterate through each section

times = zeros(num_sects,1);
for i = 1:floor(num_sects);
    tic
    % Identify the start and end indices of a section
    start_ind = sect_size/2*(i-1)+1;
    end_ind = sect_size+start_ind-1;
    % Get a segment of the original noise
    % Need to include points before the start_ind of the section, in case
    % the ITD is negative at the start and requires a sampling point before
    % start_ind
%     if (start_ind-ITDind < 0),
%         bb_interp = [bb_noise(end-ITDind:end) bb_noise(start_ind:end_ind+ITDind)]',upsamp);
%     elseif  (end_ind+ITDind > sect_size),
%         bb_interp = rect_upsample([bb_noise(start_ind-ITDind:end_ind) bb_noise(1,ITDind)]',upsamp);
%     else
%         bb_interp = rect_upsample(bb_noise(start_ind-ITDind:end_ind+ITDind)',upsamp);
%     end
    st_ups = (start_ind-1)*upsamp+1;
    end_ups = (end_ind+2*ITDind)*upsamp;
    bbsect = bb_ups(st_ups:end_ups);
    upsind = ITDind*upsamp; % index shift amount, to account for added signal of duration ITDind
    for n = 1:sect_size/2,
        try
%         CO_chan(start_ind+n-1) = bbsect((start_ind-1+n-1+ITDind)*upsamp+1);
        CO_chan(start_ind+n-1) = bbsect(upsind + (n-1)*upsamp+1);
        catch err
            error(err);
        end
        % Pick the delayed point from the upsampled signal:
        % - upsind (=ITDind*upsamp) -- number of samples in the upsampled
        % segment, equal to the maximal amount of delay in the ITD
        % trajectory
        % - (n-1)*upsamp+1 -- The sample point in the upsampled signal equal to
        % the sample point in the original noise at n
        % - ITD(start_ind+n-1,r)/10^6*Fs*upsamp -- the number of samples to
        % move in the upsampled signal, either forward or backward, in
        % order to retrieve the delayed point equivalent to the ITD. Round
        % to the closest sample
        IP_chan(start_ind+n-1) = bbsect(upsind + (n-1)*upsamp+1 + round(ITD(start_ind+n-1)/10^6*Fs*upsamp));
    end
    times(i) = toc;
end

disp(['Average elapsed time: ' num2str(mean(times)) ' +/- ' num2str(std(times)) ' s']);

%% Mix with uncorrelated noise, based on ITDdepth
uncorr = randn(2,length(t));
mix_rat = 1/(1+sqrt(1/(10^(ITDdepth/20))-1)); % compute mixing ratio
CO_chan = CO_chan*mix_rat + uncorr(1,:)*(1-mix_rat);
IP_chan = IP_chan*mix_rat + uncorr(2,:)*(1-mix_rat);

%% Apply amplitude modulation
%         env = (1 + 10^(SAMdepth/10)*cos(2*pi*t*SAMmod + phase*(2*pi)))/2;
if strcmp(AM_type,'sinusoidal'),
    env = (1 + 10^(SAMdepth/20)*cos(2*pi*t*SAMmod + 0.5*(2*pi)))/2;
elseif strcmp(AM_type,'triangle'),
    env = (1 + 10^(SAMdepth/20)*sawtooth(2*pi*t*SAMmod,width))/2;
end
CO_chan = CO_chan.*env;
IP_chan = IP_chan.*env;       

%% Save the sounds
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
%         L_file = ['DynITD\stimuli\DynITD_L' ...
%             '_-c' num2str(cntr) ...
%             '_-r' num2str(r) ...
%             '_-tok' num2str(tok)];
%         CO_file = ['DynITD\fast_stim\DynITD_CO' ...
%             '_-v' num2str(r) ...
%             '_-c' num2str(cntr) ...
%             '_-Rep' num2str(tok)];
%         IP_file = ['DynITD\fast_stim\DynITD_IP' ...
%             '_-v' num2str(r) ...
%             '_-c' num2str(cntr) ...
%             '_-Rep' num2str(tok)];
%         EPLwavwrite(CO_EPLchan,Fs,16,CO_file,'rms_dB',db(sqrt(mean(CO_EPLchan.^2))));
%         EPLwavwrite(IP_EPLchan,Fs,16,IP_file,'rms_dB',db(sqrt(mean(IP_EPLchan.^2))));
%         soundsc([CO_EPLchan; IP_EPLchan]', Fs);
%     sound(zeros(1,1:length(2*Fs)))
%         files = [files; {CO_file IP_file}];
files = [IP_EPLchan; CO_EPLchan]';