function data_clean = clean_one_channel_redundancy(data, nchan_f, npts)
%
%   function data_clean = CleanOneChannelRedundancy(data, nchan_f, npts)
%
% Description:
%   In the ONE channel case, each trial contains TWO concatenated SAME vectors & one all-zero vector.
%

% Check the inputs:
assert( 3==nchan_f, '--> Error at [CleanOneChannelRedundancy.m]: This is NOT a one channel recoring!');
%assert( length(data)==size(data,1), '--> Error at [CleanOneChannelRedundancy.m]: WRONG input (data)!')
assert( isvector(data), '--> Error at [CleanOneChannelRedundancy.m]: WRONG input (data)!')
if ~iscolumn(data)
    data = data(:);
end


%%
Ldata = length(data);

% size of one chunk of trial
Ntrl = load.trialSize(nchan_f, npts);        %Ntrl = nchan_f * npts;  
Nchunks = Ldata/Ntrl;
%Lshort = Ldata/nchan_f;     % length of the data without the redundancy

assert( fix(Ldata/Ntrl) == Nchunks, '--> Error at [CleanOneChannelRedundancy.m]: WRONG input (data)!' );

data = reshape(data, Ntrl, Nchunks);

data_clean = data(1:npts,:);
data_clean = data_clean(:);
