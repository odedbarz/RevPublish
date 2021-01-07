function [H, SE, tbl_meas, pars] = calc_raw_means(Traw, fs, binwidth, duration_sec, spikechan, fun)
%
%   function [avg, SE, Tout, pars] = calc_raw_means(Traw, fs, binwidth, duration_sec, spikechan, fun)
%
% Input:
%   Traw        : (table) contains raw waveforms & data
%   fs          : (Hz), (1x1) sampling rate of the raw waveforms in Traw
%   spikechan   : (1x1) spike-channel measurement to process 
%   fun         : (str) can be MUA or LFP
%
% Output:
%   avg         : (Nx1) means of all trials; N is the number of DRR measurement cases
%   SE          : (Nx1) standard errors of the AVG signals.
%
% Description:
%   Calculate average and standard error of the mean (SE) vectors for a given spike channel.
%

if 6 > nargin || isempty(fun)
    fun = 'MUA';    % {'MUA', 'LFP'}
end

% New downsample frequency
fs_dwn = 1/(1e-3*binwidth);

outer_list = unique(Traw.outer);
assert(~isempty(outer_list), '--> [ERROR at calc_raw_means]: no such OUTER indices!');
len_outer = length(outer_list);

inner_list = unique(Traw.inner);
assert(~isempty(inner_list), '--> [ERROR at calc_raw_means]: no such INNER indices!');
len_inner = length(inner_list);

idx.spikechan = Traw.spikechan == spikechan;
assert( 0 < nnz(idx.spikechan), ...
    '--> [ERROR at calc_raw_means]: couldn''t find the requested SPIKECHAN');

avg = cell(1, len_outer * len_inner);
SEc = cell(1, len_outer * len_inner);

Dist        = nan(len_outer*len_inner, 1);
Reverb      = nan(len_outer*len_inner, 1);
single_meas = cell(len_outer*len_inner, 1);
tbl_meas    = table(Dist, Reverb, single_meas);
clear Dist Reverb

% # of samples for a FULL measurement
duration_smp = round(duration_sec * fs);

for n = 1:len_outer
    idx.outer = outer_list(n) == Traw.outer;
    
    for m = 1:len_inner
        idx.inner = inner_list(m) == Traw.inner;
        
        % Get all trials for the current INNER & OUTER & SPIKECHAN 
        trials_nm = idx.inner & idx.outer & idx.spikechan;     % requested         
        if 0 == nnz(trials_nm)
            continue;
        end

        % Some measurements can have missing trials
        % Option #1: pad with zeros to the maximum length;
        % Option #2: throw away this measurement
        duration_x = arrayfun(@(I) length(Traw.x{I}), find(trials_nm));
        if ~all(duration_smp == duration_x)
            continue;
        end
        
        Imn = m + (n-1)*len_inner;
        [X, pars] = feval(fun, [Traw.x{trials_nm}], fs, fs_dwn);    % fun: MUA or SU
        avg{Imn} = median(X, 2);
        SEc{Imn} = median( abs(X - avg{Imn}), 2);
        
        tbl_meas.Dist(Imn)  = unique(Traw.Dist(trials_nm));
        tbl_meas.Reverb(Imn)= unique(Traw.Reverb(trials_nm));
        tbl_meas.single_meas{Imn}= X;
    end
    
end

% Transform into one matrix, and make sure to assign the vectors into the right
% columns
midx = Impale_entries_to_indices(tbl_meas);     % measurement's indices
H = nan(length(avg{1}), length(midx));
H(:,midx)  = [avg{:}];     % mean response   

SE = nan(length(avg{1}), length(midx));
SE(:,midx)  = [SEc{:}];     % STD of the mean response 







