function T = RIR(num_taps, Fs, walls, wtypes, src_angle, src_dist, sph_cent, recs_dist)
% 
%   function [Yrev, rir] = RIR(yv, Fs, walls, wtypes, src_angle, src_dist, sph_cent, recs_dist)
% 
%

assert(isscalar(num_taps) || isempty(num_taps), '--> ERROR at [RIR.m]: num_taps must be scalar!!');
if isempty(num_taps)
    % Number of samples in the RIR
    num_taps = 50e3;    % default value
end

sph_r       = 0;                    % no sphere
c_snd       = 344.5;                % (m/s) sound velocity in air   
highpass    = 0;                    % don't filter out DC component
dsply       = 1;                    % show progress
err         = 0;                    % no error introduced
Nwtypes     = length(wtypes);
Ndist       = length(src_dist);
Nangle      = length(src_angle);


fprintf('--> num_taps : %g (samples)\n', num_taps);
fprintf('--> RIR      : %g (sec)\n',     num_taps/Fs);
fprintf('--> wtypes   : %s (degrees)\n', num2str(wtypes, '%.2f '));
fprintf('--> src_angle: %s (degrees)\n', num2str(src_angle, '%.2f '));
fprintf('--> src_dist : %s (degrees)\n', num2str(src_dist, '%.2f '));

assert( isvector(wtypes) );
assert( isvector(src_angle) );
assert( isvector(src_dist) );

if ~isrow(wtypes), wtypes = wtypes'; end
if ~isrow(src_angle), src_angle = src_angle'; end
if ~isrow(src_dist), src_dist = src_dist'; end
    

%% Craete the rir table
Ntable    = Nwtypes*Ndist*Nangle;
var_names = {'wtypes', 'dist', 'angle', 'sources', 'rir', 'lead_zeros',...
    'drrL', 'drrR', 'drr_thr'};
T         = array2table( zeros(Ntable, length(var_names)), 'VariableNames', var_names);
T.rir     = cell(Ntable,1);
T.sources = cell(Ntable,1);     
T.drrL    = nan(Ntable,1);     
T.drrR    = nan(Ntable,1);     

% Populate the table with all the values in {'wtypes', 'dist', 'angle'}
all_comb = combvec(wtypes, src_angle, src_dist)';
for qq = 1:Ntable
    T.wtypes(qq) = all_comb(qq,1); 
    T.angle(qq)  = all_comb(qq,2); 
    T.dist(qq)   = all_comb(qq,3); 
end

% Add constant parameters 
T.Fs       = Fs*ones(Ntable, 1);
T.recs     = cell(Ntable, 1);   
T.walls    = cell(Ntable, 1);  T.walls(:) = {walls};
T.sph_cent = cell(Ntable, 1);  T.sph_cent(:) = {sph_cent};
T.sph_r    = sph_r*ones(Ntable, 1);   
T.c_snd    = c_snd*ones(Ntable, 1);   
T.num_taps = num_taps*ones(Ntable, 1);   
T.highpass = highpass*ones(Ntable, 1);   
T.dsply    = dsply*ones(Ntable, 1);   
T.err      = err*ones(Ntable, 1);   


%% The head 
% Two receivers 11 cm apart (11 cm rabbit; 21 cm for humans)
%
% Using two pressure sensors with no head in between (so you get ITDs but no ILDs), 
% you need an inter-sensor distance of 10.29 cm, assuming a speed of sound of 343 m/s.
recs(1,:) = sph_cent + [-recs_dist/2 0 0];
recs(2,:) = sph_cent + [recs_dist/2 0 0];

T.recs(:) = {recs};

%% Sound source

% Initialize the sources into the table
for qq = 1:Ntable
    [xt, yt, zt] = polar2rect(T.dist(qq), T.angle(qq));
    T.sources{qq,:} = sph_cent + [xt yt zt];
end


%% The RIR + convolution with the signal
for ii = 1:Ntable
    aux.cprintf('cyan', '--> (%d) Source distance : %g cm\n', ii, 100*T.dist(ii));
    aux.cprintf('cyan', '--> (%d) Source distance : %g cm\n', ii, 100*T.angle(ii));
    aux.cprintf('cyan', '--> (%d) Wall coefficient: %g cm\n', ii, T.wtypes(ii));
    
    % rir
    [h_rir, T.lead_zeros(ii)] = room_impulse(...
        T.sources{ii,:},...
        recs,...
        walls,...
        T.wtypes(ii),...
        sph_cent,...
        sph_r,...
        Fs,...
        c_snd,...
        num_taps,...
        highpass,...
        dsply,...
        err ...
    );
        
    % Normalize response to 1
    assert( 0 ~= prod(max(h_rir)), '--> There is something wrong with the RIR response!' )
    T.rir{ii} = h_rir./max(h_rir);
    
    % Direct-to-reverberation ratio
    [drr, T.drr_direct_zone(ii)] = Direct_to_Reverberation( T.rir{ii} );
    T.drrL(ii) = drr(1);
    T.drrR(ii) = drr(2);
    
end













