function tl = spike_latencies(t, ch, lmin, order)
% LATENCIES - Compute Nth order spike latencies to each sync pulse
% LATENCIES(T, CH) returns the first-order latencies 
% LATENCIES(T, CH, LMIN) ignores trials for which latency < LMIN 
% LATENCIES(T, CH, LMIN, ORDER) returns latencies of specified ORDER
%
if nargin < 3, lmin = 0; end
if nargin < 4, order = 1; end
syncchan = 0;
spikechan = 1;

% Trial markers:
ksync = find(ch == syncchan);
tl = Inf(size(ksync));

if ch(end) ~= syncchan,
   ksync(end+1) = length(ch)+1;
end


for k = 1:length(ksync)-1,
   ksp = ksync(k)+1 : ksync(k+1)-1;
   tsp = t(ksp);
   tsp = tsp(ch(ksp) == spikechan) - t(ksync(k));

   if ~isempty(tsp) && length(tsp)>=order,
      tl(k) = tsp(order);
   end
end

tl(tl<lmin) = NaN;
