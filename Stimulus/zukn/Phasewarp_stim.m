function files = Phasewarp_stim(mod_rates, num_tokens, varargin)
% Create phase warp stimuli, based on Siveke et al (2007, 2008)
% Inputs:
%   - mod_rate = modulation frequency (Hz)
%   - num_tokens = the number of tokens for which to create the stimuli
% Outputs:
%   - files = a cell array of filenames for the stimuli created

%%% (8/5/2014) Make sure that the rms value is 1 (should be for these signals,
%%% just like bb noise)

% Variables
Fs = 100000;
nbits = 16;
bandlimit = 0;
band_cent = 500;
pass_band = 400;
stim_dur = 5; % stimulus duration, in s
start_tok = randi(100);
shift_dir = 1; % specifies the direction of the phase shift
path = 'V:\ITD\Phasewarp\';
save_flag = 0;
play_flag = 0;

% Parse varargin
if ~isempty(varargin),
    for n=2:2:length(varargin),
       eval([varargin{n-1} '=varargin{n};']); 
    end
end

% Create stimuli
files = {}; % to save filenames
for m = 1:length(mod_rates),
    for n = 1:num_tokens,
        rng(start_tok+n-1);
        if bandlimit,
            % If desired, zero out frequencies outside the passband
            pass_inds = pass_band*stim_dur; % number of indexes in the pass band
            cent_ind = band_cent*stim_dur;
            mag = zeros(Fs*stim_dur/2,1);
            mag(cent_ind-round(pass_inds/2):cent_ind+round(pass_inds/2)) = 1;
            mag = [mag; flipud(mag)];
        else
            % Create noise with constant magnitude and random phase
            mag = ones(Fs*stim_dur,1);
        end
        
        ph_fst = rand(Fs*stim_dur/2,1)*2*pi; % Phase in first half of spectrum
        ph_shf_fst = circshift(ph_fst,shift_dir*mod_rates(m)*stim_dur); 
            % The phase is the first half of the spectrum is circularly shifted by
            % the modulation frequency
        ph_snd = flipud(ph_fst); % Phase in the second half of the spectrum is a flipped
            % version of the first half (in order to make the spectrum an odd spectrum,
            % true for real-valued signals)
        ph_shf_snd = flipud(ph_shf_fst);
        phase = [ph_fst; ph_snd]; % Non-shifted phase array
        phase_shf = [ph_shf_fst; ph_shf_snd]; % Shifted phase array
        
        L = ifft(mag.*exp(1i*phase));
        R = ifft(mag.*exp(1i*phase_shf));
        y = real([L, R]); % Take the real component
        % Save stimuli
        mu = 1/max(max(abs(y))); % to normalize y by the stimuli
        y = y.*mu;
        if save_flag,
            IP_file = [path 'Phasewarp_IP_-v' num2str(mod_rates(m)) '_-Rep' num2str(n)];
            CO_file = [path 'Phasewarp_CO_-v' num2str(mod_rates(m)) '_-Rep' num2str(n)];
            EPLwavwrite(y(:,1),Fs,nbits,IP_file,'rms_dB',db(sqrt(mean(y(:,1).^2))));
            disp(IP_file);
            EPLwavwrite(y(:,2),Fs,nbits,CO_file,'rms_dB',db(sqrt(mean(y(:,2).^2))));
            disp(CO_file);
            files = [files {IP_file, CO_file}];
        end
        if play_flag,
            soundsc(y,Fs);
        end
        files = y;
    end
end

%soundsc(y,Fs); % play stim