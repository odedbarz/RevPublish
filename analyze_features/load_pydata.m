function [pyspec, pypitch, pydata] = load_pydata(n_time)
% pydata = 
%   struct with fields:
% 
%     stim: [1×1 struct]
%        p: [1×1 struct]
%     harm: [1×1 struct]
%     spec: [1×1 struct]
mat_fn      = 'analyzed_librosa.mat';
mat_full_fn = fullfile( load.path_to_data('data'), 'Analysis', mat_fn);
pydata      = load( mat_full_fn );

% Extract parameters
% pyspec = pydata.spec;
pypitch      = pydata.p;

pyspec      = pydata.spec;
pyspec.f0_smp = 5;  % don't start from f == 0 Hz
pyspec.f    = pyspec.f(pyspec.f0_smp:end);     % Remove low frequencies
pyspec.Sdb  = pyspec.Sdb(pyspec.f0_smp:end,:);
pyspec.S    = pyspec.S(pyspec.f0_smp:end,:);

% Interpulate to N_TIME
pypitch.F0i = interp1( linspace(1,pypitch.duration_seconds,pypitch.nt_), pypitch.F0,...
    linspace(1,pypitch.duration_seconds,n_time) )';
pyspec.Sdbi = interp1( linspace(1,pypitch.duration_seconds,pypitch.nt_), pyspec.Sdb',...
    linspace(1,pypitch.duration_seconds,n_time) )';
pyspec.Si   = interp1( linspace(1,pypitch.duration_seconds,pypitch.nt_), pyspec.S',...
    linspace(1,pypitch.duration_seconds,n_time) )';
pyspec.ti   = interp1( linspace(1,pypitch.duration_seconds,pypitch.nt_), pyspec.t,...
    linspace(1,pypitch.duration_seconds,n_time) )';
