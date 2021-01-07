% Let i be the vector of spike times
% Let j be the vector of trial numbers
% Let T be the length of each trial (in same time units as times in <i>)
% Let <BinSize> be a user specified scalar defining desired bin size (same
%     units as spike times in <i>)

% Define a new vector <SpTimes> where the spike times are
% redefined relative to the beginning of trial 1 (as opposed
% to being defined relative to the beginning of each trial) giving
% time stamp values 'Post-start-time' rather than 'post-stimulus-time'
SpTimes = i + (j-1)*T

% First define number of bins for histogram to acheive desired bin size
NumberofBins = max(SpTimes)./BinSize;

% Now use Matlab HIST function to find histogram and bin numbers
[NumberofSpikes BinPosition] = hist(SpTimes,NumberofBins)

% Find total length number of bins from beginning of first trial to end
% of last trial
TotalBins = ceil(T./BinSize)*max(j);

% Create a vector of length <TotalBins> full of zeros
PostStartTimeHistogram = zeros(1,TotalBins)

% For each value of <BinPosition>, put corresponding value from <NumberofSpikes>
PostStartTimeHistogram(BinPosition) = NumberofSpikes;

% Now reshape matrix to get raster
Raster = reshape(PostStartTimeHistogram,ceil(T./BinSize),max(j));

