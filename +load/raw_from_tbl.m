function [x, t] = raw_from_tbl(raw_st, outer, inner, channel, trial, verbose)

if 6 > nargin
    verbose = false;
end
%channel   = 3;    %1+fix(spikechan/6);     % channels number
%spikechan = 1 + (channel-1)*6;

% Find the desired row in the table
idx.channel = raw_st.tbl.channel == channel;
idx.trial   = raw_st.tbl.trial == trial;
idx.inner   = raw_st.tbl.inner == inner;
idx.outer   = raw_st.tbl.outer == outer;
idx.meas    = idx.channel & idx.trial & idx.inner & idx.outer;
assert(1==nnz(idx.meas), '--> Can''t find the required measurement!');
% meas_row = find(idx.meas);

if verbose
    fprintf('--> # channel: %d\n', nnz(idx.channel));
    fprintf('--> # trial  : %d\n', nnz(idx.trial));
    fprintf('--> # inner  : %d\n', nnz(idx.inner));
    fprintf('--> # outer  : %d\n', nnz(idx.outer));
    fprintf('--> # meas  : %d\n', nnz(idx.meas));
end

x = raw_st.tbl.x{idx.meas};


duration_ms = raw_st.rawdata.t_duration;
duration_sec = 1e-3*duration_ms;
len_x = raw_st.sr * duration_sec;
assert( len_x == length(x) );
t = linspace(0, duration_sec, length(x))';   % (sec) time axis
