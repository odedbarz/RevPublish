function x = raw_extract_meas(X, rawdata, varargin)
%
%   function x = load.raw_extract_meas(X, rawdata, [LMA_channel], [trial],...
%       [inner], [outer], [fignum])
%


p = inputParser;

% === Required inputs ===
addRequired(p, 'X', @isnumeric);            % (smp x 1) raw vector
addRequired(p, 'rawdata', @isstruct);       % raw data structure

% === Optional inputs ===
addOptional(p, 'LMA_channel', 1, @isnumeric);  % LMA_channel = 1+fix(spikechan/6);           
addOptional(p, 'spikechan', [], @isnumeric);	% spikechan = 1 + (LMA_channel-1)*6;         
addOptional(p, 'trial', 1, @isnumeric);           
addOptional(p, 'inner', 1, @isnumeric);           
addOptional(p, 'outer', 1, @isnumeric);           
addOptional(p, 'fignum', [], @isnumeric);           

parse(p, X, rawdata, varargin{:});

LMA_channel = p.Results.LMA_channel;       
spikechan   = p.Results.spikechan;       
trial       = p.Results.trial;
inner       = p.Results.inner;       
outer       = p.Results.outer;
fignum      = p.Results.fignum;       

% In this case, over-write the LMA_channel
if any(contains(p.UsingDefaults, 'LMA_channel')) && ~isempty(spikechan) 
    LMA_channel = 1 + fix(spikechan/6);  
end

assert(rawdata.n_spikechan >= LMA_channel, '--> [+load.raw_extract_meas]: wrong LMA_channel number!');


%%
%N = rawdata.smp_token;
n_spikechan = rawdata.n_spikechan;
smp_token   = rawdata.smp_token;
n_tokens    = rawdata.n_tokens;
n_trials_per_stimulus = rawdata.n_trials_per_stimulus;
n_inner     = max(1,length(rawdata.innerSeq));



%% Impale measurement ordering
% Get the right inner & outer measurements; use Impale recording order.
inner = find(rawdata.inner_iorder == inner);
assert(isscalar(inner));

outer = find(rawdata.outer_iorder == outer);
assert(1 == length(outer));



%%
% Reshape the RAW vector into tokens (columns)
X_ = reshape(X, smp_token, []);  

bias.spikechan = n_spikechan*[(1:n_tokens)-1];
bias.channel   = (LMA_channel-1);
bias.trial     = (trial-1)*n_tokens*n_spikechan;
bias.inner     = (inner-1)*n_trials_per_stimulus*n_tokens*n_spikechan;
bias.outer     = (outer-1)*n_inner*n_trials_per_stimulus*n_tokens*n_spikechan;
meas_cols      = 1 + bias.outer +  + bias.inner + bias.trial + bias.channel + bias.spikechan;

x = X_(:, meas_cols);
x = x(:);


%% Polarity
x = rawdata.polarity(LMA_channel) * x;



%% Plot
if isempty(fignum)
    return;
end

fprintf('--> rawdata:\n');
disp(rawdata);

figure(99);
clf;
%Fs = 100e3;     % (Hz) Impale's default sampling rate
t = linspace(0, 1e-3*rawdata.t_duration, rawdata.smp_duration)';  
plot(t, [x, abs([0; diff(x)])]);
xlabel('Time (sec)');
legend('x', 'diff(x)');




