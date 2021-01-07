function Ntrl = trialSize(nchan_f, npts)
%
% Number of samples in each trial in the data vector. For one channel, each
% trial contains two concatenated signals and one all-zero vector.
%

Ntrl = nchan_f * npts;  % samples