function [files] = MotSource_ITDnAM(cntr_ITD, amp_ITD, mod_rates, SAMmod, SAMdepth, ITDdepth, phase, varargin)
%%% Creates stimuli with both a time-varying ITD and sinusoidally varying
%%% amplitude.  The creation of the moving ITD is identical to
%%% MotSource_ITD.m.  The sinusoidal envelope is applied afterwards by
%%% multiplying the L and R channels by the envelope.  The stimulus is
%%% played with soundsc before it is saved.
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

%%% EDIT (7/28/14): Should alternate between ipsi and contra stimuli for varying ITD.
%%% This should be set so the delay is always decreasing for the side with
%%% the varying delay. (+ ITD is toward contra, - ITD = away from contra)
%%% EDIT (7/28/14): IP is in first column, CO is in second.  For
%%% psychophysics, a positive ITD is to the right, negative is to the left.
%%% This way, when the stimuli are played directly to the listener, the
%%% positive and negative ITDs are presented with the proper orientation.

%% Initial values
Fs = 100000; % sampling rate (Hz)
c = 340.29; % speed of sound in m/s (taken from google search)
% tok = 42;
tok = randi(100);

t_end = 1; % length of time for the stimulus (s)

files = {};
width = 0.5;

step_ITD = round(1/Fs*10^6); % ITD step size during ITD motion (us)
motion_type = 'triangle';
AM_type = 'sinusoidal';

if ~isempty(varargin),
    for n = 2:2:length(varargin),
        eval([varargin{n-1} '=varargin{n};']);
    end
end

for r = mod_rates
    %% Calculate the trajectory of motion
    %motion = sawtooth(2*pi*r*t,0.5);
    t = 0:1/Fs:t_end-1/Fs;
    if r,
        if strcmp(motion_type,'triangle')
            motion = sawtooth(2*pi*r*t + (0.5+phase)*(2*pi),width);
        elseif strcmp(motion_type,'sinusoidal')
            motion = cos(2*pi*r*t + phase*(2*pi));
        end
    else
        motion = zeros(1,length(t));
    end
    for cntr = cntr_ITD
%         ITD = amp_ITD*motion + cntr_ITD; % convert to ITDs
        ITD = -amp_ITD*motion + cntr;
        
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
            CO_chan = zeros(1,length(bb_noise)); %%% EDITED: 7/26/14
            IP_chan = zeros(1,length(bb_noise));

            % Iterate through ITD steps, designated by dITD_ind
            for i = 1:length(dITD_ind)-1;
                start_ind = dITD_ind(i)+1;
                end_ind = dITD_ind(i+1);
                ITD_i = ITD_st(start_ind); % get the ITD at the index
                    % corresponding to t_incr(i)
                ITD_sf = round(ITD_i/10^6*Fs); % compute the number of indeces for the shift
%                 ITD_sf = ITD_i/10^6*Fs; 
                if (start_ind-abs(ITD_sf) < 1) || (end_ind+abs(ITD_sf) > length(bb_noise)),
                    bb_wrapped = [bb_noise(length(bb_noise)-abs(ITD_sf)+1:end) ...
                        bb_noise bb_noise(1:abs(ITD_sf))];
                    st_i = start_ind+abs(ITD_sf); % correct for wrapping
                    end_i = end_ind+abs(ITD_sf);
                else
                    bb_wrapped = bb_noise;
                    st_i = start_ind;
                    end_i = end_ind;
                end
%                 if dITD(dITD_ind(i)) > 0, %%% EDITED 9/30/14
                if dITD(dITD_ind(i)) < 0, % contra experiencing decreasing delay
                    IP_chan(start_ind:end_ind) = bb_wrapped(st_i-ITD_sf:end_i-ITD_sf);
                    CO_chan(start_ind:end_ind) = bb_wrapped(st_i:end_i);
                else % ipsi experiencing decreasing delay
                    CO_chan(start_ind:end_ind) = bb_wrapped(st_i-(-ITD_sf):end_i-(-ITD_sf));
                    IP_chan(start_ind:end_ind) = bb_wrapped(st_i:end_i);
                end
%                 IP_chan(start_ind:end_ind) = bb_wrapped(st_i-ITD_sf:end_i-ITD_sf);
            end

        %% Mix with uncorrelated noise, based on ITDdepth
        uncorr = randn(2,t_end*Fs);
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
    end
end
