function [anf, cfs, t, Fs] = Stim2ANF( stim, varargin)
%
% function [anf, cfs, t, Fs] = Stim2ANF( stim, [Nch], [downsmp], [lowfreq], [highfreq])
%
%   stim: (Tx1) the stimulus to convert to auditory nerve response (ANF);
%         !!! NOTE: <stim> must be sampled at 100k Hz
%
%
% Description:
%   Converts the stimulus stim to ANF activity.
%
% Example:
%
%   Run Meddis' model with x1 tone od 2.5k Hz:
%   >> Stim2ANF( 0.001*sin(2*pi*[0:0.1*100e3-1]/100e3 * 2500), 'Meddis', 'fignum', 11);
%
%   Run Meddis' model with Fs = 25k Hz:
%
%       (Fs = 100kHz)
%       sfun = @(ff) ones(1,length(ff)) * 0.001*sin(2*pi*[0:0.1*100e3-1]/100e3 .* ff(:));
%       [anf, cfs, Fs] = Stim2ANF( sfun(1500) + sfun(2750), 'fignum', 11);
%
%   or for the Zilany & Carney's model,
%
%   	[anf, cfs, Fs] = Stim2ANF( sfun(1500)+sfun(2750), 'method', 'Carney', 'fignum', 11);
%
%

%% Parse the input
p = inputParser;

addRequired(p, 'stim', @isnumeric);                % (1xT) the stimulus to convert to anf

chl_method = @(method)...
    any(cellfun(@(SS) strcmpi(SS,method), {'none', 'carney', 'slaney', 'meddis', 'gammatone'}, 'UniformOutput', 1 ));
addOptional(p, 'method', 'meddis', chl_method);           % (str) the method to use {'Zilany', 'Meddis'}

addOptional(p, 'Fs', 100e3, @isnumeric);        % [Hz] stimulus' sampling rate
addOptional(p, 'Nch', 256, @isnumeric);         % # of channels\filters in cochlea
addOptional(p, 'downsmp', 0, @isnumeric);       % downsampling factor
addOptional(p, 'lowfreq', 360.0, @isnumeric); 	% [Hz] lowest frequency (for rabbits)
addOptional(p, 'highfreq', 42e3, @isnumeric); 	% [Hz] highest frequency (for rabbits)
addOptional(p, 'align_cfs', 'ascend', @(s) any(strcmpi({'ascend', 'descend'}, s)) ); % arrange the CFs

% gammatone
addOptional(p, 'calc_envelope', false, @(x) islogical(x) || isscalar(x));


% Meddis' parameters:
addOptional(p, 'rabbit_flag', false, @islogical);% use rabbit's cochlea parameters

% Zilany's parameters:
addOptional(p, 'cohc', 1.0, @isnumeric);        % Condition (normal->bad) of ohc function
addOptional(p, 'cihc', 1.0, @isnumeric);        % Condition (normal->bad) of ihc function
addOptional(p, 'species', 1, @isnumeric);       % 1 for cat (2 for human with Shera et al. tuning; 3 for human with Glasberg & Moore tuning)
addOptional(p, 'noiseType', 0, @isnumeric);     % 1 for variable fGn (0 for fixed fGn)
addOptional(p, 'fiberType', 3, @isnumeric);     % spontaneous rate (in spikes/s) of the fiber BEFORE refractory effects; "1" = Low; "2" = Medium; "3" = High
addOptional(p, 'implnt', 0, @isnumeric);        % "0" for approximate or "1" for actual implementation of the power-law functions in the Synapse
addOptional(p, 'nrep', 1, @isnumeric);          % The number of stimulus repetitions, e.g., 50 (always 1 by default in this implementation).

addOptional(p, 'fignum', [], @isnumeric);       % figure number to plot (if so desired)

parse(p, stim, varargin{:});

stim        = p.Results.stim;       % (--> row vector)
method      = p.Results.method;
align_cfs   = p.Results.align_cfs;
calc_envelope= p.Results.calc_envelope;
Fs          = p.Results.Fs;
rabbit_flag = p.Results.rabbit_flag;
Nch         = p.Results.Nch;
downsmp     = p.Results.downsmp;
lowfreq     = p.Results.lowfreq;
highfreq    = min(Fs, p.Results.highfreq);
fignum      = p.Results.fignum;

% Make sure that <stim> is a vector
assert( iscolumn(stim) || isrow(stim),...
    '--> [ERROR in Stim2ANF]: the input stimulus <stim> MUST be a vector!')

% Make sure that <stim> is a row vector
if iscolumn(stim), stim = stim'; end

% The CF filters in the cochlea
cfs = ERBSpace(lowfreq, highfreq, Nch, rabbit_flag);     % (also implemented in MakeERBFilters.m)



%%
switch lower(method)
    case 'none'
        % Stimuli WITHOUT periphery
        cfs = ERBSpace(lowfreq, highfreq, Nch, rabbit_flag);

        Nt_ = length(stim);
        stimulus_time = Nt_*1/Fs;   % (sec) total stimulus time
        t_ = [1/Fs:1/Fs:stimulus_time];

        stimulus_fun = @(FF) ( ones(length(FF),1) * sin(2*pi*t_'.*FF(:)) )';

        anf = cell2mat(arrayfun(@(X) stimulus_fun(X), cfs, 'UniformOutput', false) );

    case 'gammatone'
        % Slaney's toolbox (see AuditoryToolboxTechReport)
        fcoefs = MakeERBFilters(...
            Fs,...          % Sampling rate
            Nch,...         % number of channels in the cochlea
            lowfreq, ...  	% lowest frequency in the cochlea to simulate
            highfreq, ...  	% highest frequency in the cochlea to simulate
            rabbit_flag ... % use rabbit's cochlea parameters
        );
        coch = ERBFilterBank( stim, fcoefs );

        if calc_envelope
            anf = envelope(coch')';
        else
            anf = coch;
        end

    case 'slaney'
        % Slaney's toolbox (see AuditoryToolboxTechReport)
        % MIDDLE EAR filter\Pre-emphasis filter (Meddis & O'Mard, 1997)
        ear_filter = OuterEarFilt(Fs);
        s_filt = filter( ear_filter, stim );

        fcoefs = MakeERBFilters(...
            Fs,...          % Sampling rate
            Nch,...         % number of channels in the cochlea
            lowfreq, ...  	% lowest frequency in the cochlea to simulate
            highfreq, ...  	% highest frequency in the cochlea to simulate
            rabbit_flag ... % use rabbit's cochlea parameters
        );
        coch = ERBFilterBank( s_filt, fcoefs );

        % From the AuditoryToolboxTechReport, page 24-25:
        % Note: "...there is no adaptation or automatic gain control to
        %       equalize the formant frequencies and enhance the onsets."
        %
        % (below is my implementation using MATLAB's vectorization)
        c = max(coch, 0);                   	% half-wave rectifier
        anf = filter(1, [1 -.99], c, [], 2); 	% low-pass filter


    case 'meddis' 	% Ray Meddis� 1986 JASA paper hair cell model
                    % (see also AuditoryToolboxTechReport)
        % MIDDLE EAR filter\Pre-emphasis filter (Meddis & O'Mard, 1997)
        ear_filter = OuterEarFilt(Fs);
        s_filt = filter( ear_filter, stim );

        fcoefs = MakeERBFilters(...
            Fs,...          % Sampling rate
            Nch,...         % number of channels in the cochlea
            lowfreq, ...  	% lowest frequency in the cochlea to simulate
            highfreq, ...  	% highest frequency in the cochlea to simulate
            rabbit_flag ... % use rabbit's cochlea parameters
        );
        coch = ERBFilterBank( s_filt, fcoefs );
        anf = MeddisHairCell( 80* coch/max(coch(:)), Fs);


    case 'carney'  	% see Zilany & Carney's (2014)Code_and_paper file
        % Carney's model parameters:
        cohc        = p.Results.cohc;       % ==1 for normal ohc function
        cihc        = p.Results.cihc;       % ==1 for normal ihc function
        species    	= p.Results.species;    % 1 for cat (2 for human with Shera et al. tuning; 3 for human with Glasberg & Moore tuning)
        noiseType  	= p.Results.noiseType;  % 1 for variable fGn (0 for fixed fGn)
        fiberType  	= p.Results.fiberType;  % spontaneous rate (in spikes/s) of the fiber BEFORE refractory effects; "1" = Low; "2" = Medium; "3" = High
        implnt      = p.Results.implnt;     % "0" for approximate or "1" for actual implementation of the power-law functions in the Synapse
        nrep        = p.Results.nrep;
        %psthbinwidth = params.anf.zilany.psthbinwidth; % binwidth in seconds;

        % If the <rabbit_flag> is ON, then teh cat's model should be used
        % (species = 1); the cat's cochlea fits better to the rabbit's cochlea!
        if rabbit_flag && 1 ~= species
            warning('--> [Warning at Stim2ANF:] <rabbit_flag> is ON; Setting species to 1 (cat model)!!!');
            species = 1;
        end

        % Note: For Carney's model:
        %   * CF must be between 125 Hz and 20 kHz for HUMAN model
        %   * CF must be between 125 Hz and 40 kHz for CAT model
        assert(lowfreq >= 125, 'lowfreq MUST be bigger or equal than 125!');
        %lowfreq     = max(125, lowfreq);    % (Hz)
        %highfreq    = min(40e3, highfreq);  % (Hz)

        %cfs = ERBSpace( lowfreq, highfreq, Nch); 	% cochlear frequencies

        assert(100e3 == Fs, '--> ERROR in [Stim2ANF.m]: For the Zilany & Carney model, Fs MUST be 100k Hz!');
        T  = length(stim)*(1/Fs);  % total stimulus' time
        Ts = 1/Fs;
        Ls = length(stim);         % (samples) stimulus' length
        anf = zeros(Nch, Ls);

        warning('off', 'MATLAB:RandStream:ActivatingLegacyGenerators');
        warning('off', 'MATLAB:RandStream:ReadingInactiveLegacyGeneratorState') ;

        parfor ch = 1:Nch
            CF = cfs(ch);

            % calc. the IHCs voltage at each CF:
            % * (secs+Ts): avoid the numerical issue of (reptime < pxbins*tdres) over (reptime <= pxbins*tdres).
            %   make sure that <reptime> is greater than the stimulus duration
            % * species:
            %   (1) for cat;
            %   (2) for human with Shera et al. tuning;
            %   (3) for human with Glasberg & Moore tuning.
            dummy = model_IHC( stim, CF, nrep, Ts, T+2*Ts, cohc, cihc, species );
            vihc = dummy(1:Ls);   % remove the last added sample.

            % from IHCs voltage to mean rate:
            [meanrate, ~, ~] = model_Synapse( vihc, CF, nrep, Ts, fiberType, noiseType, implnt );
            %[meanrate, varrate, psth] = model_Synapse( vihc, CF, nrep, Ts, fiberType, noiseType, implnt );

            % Use this case for IHC fibers with various rates
            %{
            [highsr,~,~]    = model_Synapse( vihc, CF, nrep, Ts, 3, noiseType, implnt );    % fiberType == 3
            [mediumsr,~,~]  = model_Synapse( vihc, CF, nrep, Ts, 2, noiseType, implnt );    % fiberType == 2
            [lowsr,~,~]     = model_Synapse( vihc, CF, nrep, Ts, 1, noiseType, implnt );    % fiberType == 1
            meanrate        = 0.6*highsr + 0.25*mediumsr + 0.15*lowsr;
            %}
            anf(ch,:) = meanrate;

        end

        % flip the CF axis (1 dim) so that the lowest CFs will be at
        % the "lower" section of the anf matrix:
        %anf = anf(end:-1:1,:);

    otherwise
        error('--> [ERROR in Stim2ANF]: Your choice of <method> is not supported!!!');
end

%% Align the CFs
% * 'reverse' is the default
if strcmpi(align_cfs, 'ascend')
    [cfs, sort_cfs_idx] = sort(cfs, 'ascend');
    anf = anf(sort_cfs_idx, :);
end


%% Downsampling
if 0 < downsmp
    Fs = Fs/downsmp;
    anf = resample(anf', 1, downsmp)';
end


%%
Nt = size(anf,2);
stimulus_time = Nt*1/Fs;   % (sec) total stimulus time
t = (1/Fs:1/Fs:stimulus_time);


%% PLOT the ANF using log scale for the frequencies:
if ~isempty(fignum)
    %stim_dwn = downsample(stim, downsmp);    % downsampled signal
    stim_dwn = resample(stim, 1, downsmp)';    % downsampled signal

    figure(fignum);
    clf;
    %set(gcf, 'Position', [196 209 1044 769])
    subplot(1,10,1:8);
    T  = length(stim_dwn)*(1/Fs);  % total stimulus' time
    ttt = 1e3* linspace(0, T-1/Fs, size(anf,2));
    fff = log(cfs);
    [TTT,FFF] = meshgrid(ttt, fff);
    surface(TTT,FFF,log(abs(anf)),'EdgeColor','none');
    ax(1) = gca;
    ax(1).YScale = 'linear';
    ax(1).YDir = 'normal';
    axis tight
    ax(1).YTickLabel = num2str(1e-3*exp(str2double(ax(1).YTickLabel)),'%.2f');
    set(gca, 'fontsize', 22);
    xlabel('Time [ms]');
    ylabel('CF [kHz]');
    title('AN Population Response');

    subplot(1,10,9:10);
    plot(sum(anf,2), fff);
    axis tight
    ax(2) = gca;
    ax(2).YScale = 'linear';
    ax(2).YDir = 'normal';
    ax(2).YTickLabel = num2str(1e-3*exp(str2double(ax(2).YTickLabel)),'%.2f');

    grid on
    axis tight
end
