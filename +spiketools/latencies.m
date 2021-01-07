function tl = latencies(t, ch, syncchan, spikechan, order)
% LATENCIES - Compute Nth order spike latencies to each sync pulse
% LATENCIES(T, CH, SYNCCHAN, SPIKECHAN) returns the first-order latencies 
% of the spikes in SPIKECHAN relative to the pulses in SYNCCHAN
%
% LATENCIES(T, CH, SYNCCHAN, SPIKECHAN, ORDER) returns the latencies of the
% specified ORDER
%
if nargin < 5, order = 1; end
if nargin < 4, spikechan = 1; end
if nargin < 3, syncchan = 0; end

ksync = find(ch == syncchan);
tsync = t(ksync);
kspike = find(ch == spikechan);
dmax = max(kspike(order+1:end)-kspike(1:end-order)) - 1;

tl = Inf;
tl = tl(ones(size(tsync)));
current_order = zeros(size(ksync));

for d = 1:dmax;
    ktest = find(ksync+d <= kspike(end) & current_order < order);
    current_order(ktest) = current_order(ktest) + (ch(ksync(ktest)+d) == spikechan);
    kgood = ktest(find(current_order(ktest) == order));
    tl(kgood) = t(ksync(kgood)+d) - tsync(kgood);
    if all(current_order(ktest) >= order), break; end
end


