%
% chk_RIRs.m
% 
% A simple function to check the RIRs.
%

clc;

% Loads TRG structure
if ~exist('trg', 'var')
    load('\\apollo\research\ENT\Delgutte\obz\!_Reverberation_Stimuli\Project 2\!trg_M.mat');
end

% Loads TRIR table
if ~exist('Trir', 'var')
    load('\\apollo\research\ENT\Delgutte\obz\2Bertrand\Project2\DRR(-12dB)_Dist(3m)\tbl_RIR_Dist(3m).mat');
end


row = 28;   % [wtype: 0.11, dist: 3.0 m; angle: 0]
% row = 27;   % [wtype: 0.44, dist: 3.0 m; angle: 0]
disp(Trir(row,:));

% A simple window to truncate the tail
n_smp = size(Trir.rir{1}, 1);
W = @(a,b) 0.5*(1-tanh(b*((1:n_smp)'-a)/n_smp));















