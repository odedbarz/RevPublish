function h = hcompute(h, t, ch, index)
% HCOMPUTE - Compute a histogram from time events
% H = HCOMPUTE(HD, T, CH, INDEX) computes a histogram using the
% parameters specified in HD and the event times/channels in T/CH
% INDEX is a vector of indices for use only with period histograms in sequences
%
if strcmp(h.type, 'none'), return; end
if length(h.value) > 0,
    error('Can''t compute neurogram from spike train');
end

if ~strcmp(h.gate.type, 'none'),
    [t, ch] = gate_events(t, ch, h.gate, h.syncchan);
end

if nargin < 4 || length(h.period) <= 1, index = 1; end

switch h.type
case 'PST'
    [tdata, nsync] = pstimes(t, ch, h.syncchan);
    kspike = find(ch == h.spikechan);
    h.data = tdata_to_1Dhist(h, tdata(kspike));
    h.synccount = nsync;
    h.spikecount = length(kspike);
case 'Interval'
    tdata = intervals(t, ch, h.spikechan, h.order);
    h.data = tdata_to_1Dhist(h, tdata);
    % ODED, 03/21/2020: if tdata is NaN (no spikes) set spikecount to zero
    %h.spikecount = length(tdata)+h.order;                              % ODED, 03/21/2020 
    h.spikecount = (length(tdata)>1) * ( length(tdata) + h.order );     % ODED, 03/21/2020
    h.synccount = sum(ch == h.syncchan);
case 'Autocorrelation'
    dim = which_dim(h);
    maxinterv = h.offset(dim) + size(h.data,dim) * h.binwidth(dim);
    tdata = all_order_interv(t, ch, h.spikechan, maxinterv);
    h.data = tdata_to_1Dhist(h, tdata);
    h.synccount = sum(ch == h.syncchan);
    h.spikecount = sum(ch == h.spikechan);
case 'Crosscorrelation'
    dim = which_dim(h);
    maxinterv = max(abs([h.offset(dim) ...
                         h.offset(dim) + size(h.data,dim)*h.binwidth(dim)]));
    tdata = cross_intervals(t, ch, h.spikechan, maxinterv);
    h.data = tdata_to_1Dhist(h, tdata);
    h.synccount = sum(ch == h.syncchan);
    h.spikecount = count_events(ch,h.spikechan);

case 'Recurrence'
    dim = which_dim(h);
    if h.offset(dim) >= 0, 
        direction = 'forward';
    elseif h.offset(dim) + size(h.data,dim)*h.binwidth(dim) <= 0, 
        direction = 'backward';
    else
        direction = 'both';
    end
    tdata = recurrence_times(t, ch, h.spikechan, direction);
    h.data = tdata_to_1Dhist(h, tdata);
    h.synccount = sum(ch == h.syncchan);
    h.spikecount = count_events(ch,h.spikechan);

case 'ShuffledCorrelation'
    tdata = shuffled_intervals(t, ch, h.syncchan, h.spikechan);
    h.data = tdata_to_1Dhist(h, tdata);
    h.synccount = count_events(ch,h.syncchan);
    h.spikecount = count_events(ch,h.spikechan);
case 'Period'
    [tdata, nperiod] = periodtimes(t, ch, h.period(index), h.syncchan);
    kspike = find(ch == h.spikechan);
    dim = which_dim(h);
%     h.binwidth(dim) = h.period(index)/length(h.data);
%     h.data = tdata_to_1Dhist(h, tdata(kspike));   
    h.binwidth(dim) = 1/length(h.data);
    h.data = tdata_to_1Dhist(h, tdata(kspike)/h.period(index));   
       
    h.synccount = sum(ch == h.syncchan); 
    h.spikecount = length(kspike);
    h.periodcount = nperiod;
case 'Latency'
    tdata = latencies(t, ch, h.syncchan, h.spikechan, h.order);
    h.data = tdata_to_1Dhist(h, tdata);   
    h.synccount = length(tdata); 
    h.spikecount = sum(ch == h.spikechan);
case 'Revcor'
    h = revcor(h, t, ch);
case 'PESE'
    h = pese(h, t, ch);
case 'PST Raster'
    [tdata, nstim] = pstimes(t, ch, h.syncchan);
    trialno = cumsum(ch == h.syncchan);
    kspike = find(ch == h.spikechan);
    h.data = tdata_to_raster(h, trialno(kspike), tdata(kspike));
    h.spikecount = length(tdata(kspike));
    h.synccount = trialno(end);
case 'Interval Raster'
    [tdata, kspike] = intervals(t, ch, h.spikechan, h.order);
    trialno = cumsum(ch == h.syncchan);
    h.data = tdata_to_raster(h, trialno(kspike(:,2)), tdata);
    h.spikecount = length(tdata)+h.order;
    h.synccount = trialno(end);   
case 'Period Raster'
    [tdata, nperiod] = periodtimes(t, ch, h.period(index), h.syncchan);
    trialno = cumsum(ch == h.syncchan);
    kspike = find(ch == h.spikechan);
    h.binwidth(2) = 1/size(h.data,2);
    h.data = tdata_to_raster(h, trialno(kspike), tdata(kspike)/h.period(index));
    h.synccount = trialno(end);   
    h.spikecount = length(kspike);
    h.periodcount = nperiod;
case 'Autocorrelation Raster'
    maxinterv = h.offset(2) + size(h.data,2) * h.binwidth(2);
    [tdata, kspike] = all_order_interv(t, ch, h.spikechan, maxinterv);
    trialno = cumsum(ch == h.syncchan);
    if ~isempty(kspike),
        h.data = tdata_to_raster(h, trialno(kspike(:,2)), tdata);
    end
    h.synccount = trialno(end);   
    h.spikecount = sum(ch == h.spikechan);
case 'Crosscorrelation Raster'
    maxinterv = max(abs([h.offset(2) ...
                         h.offset(2) + size(h.data,2)*h.binwidth(2)]));
    [tdata, kspike] = cross_intervals(t, ch, h.spikechan, maxinterv);
    trialno = cumsum(ch == h.syncchan);
    if ~isempty(kspike),
        h.data = tdata_to_raster(h, trialno(kspike(:,2)), tdata);
    end
    h.synccount = trialno(end);   
    h.spikecount = count_events(ch,h.spikechan);
  
case 'PST-Interval'
    [xdata, nsync] = pstimes(t, ch, h.syncchan);
    [ydata, kspike] = intervals(t, ch, h.spikechan, h.order);
    h.data = tdata_to_2Dhist(h, ydata, xdata(kspike(:,2)));
    h.synccount = nsync;
    h.spikecount = length(kspike) + h.order;
case 'PST-Period'
    [xdata, nsync] = pstimes(t, ch, h.syncchan);
    [ydata, nperiod] = periodtimes(t, ch, h.period(index), h.syncchan);
    kspike = find(ch == h.spikechan);
    h.binwidth(1) = 1/size(h.data,1);
    h.data = tdata_to_2Dhist(h, ydata(kspike)/h.period(index), xdata(kspike));
    h.spikecount = length(kspike);
    h.synccount = nsync;
    h.periodcount = nperiod;
case 'PST-Autocorrelation'
    [xdata, nsync] = pstimes(t, ch, h.syncchan);
    maxinterv = h.offset(1) + size(h.data,1) * h.binwidth(1);
    [ydata, kspike] = all_order_interv(t, ch, h.spikechan, maxinterv);
    if ~isempty(kspike),
        h.data = tdata_to_2Dhist(h, ydata, xdata(kspike(:,2)));
    end
    h.spikecount = sum(ch == h.spikechan);
    h.synccount = nsync;
case 'PST-Crosscorrelation'
    [xdata, nsync] = pstimes(t, ch, h.syncchan);
    maxinterv = max(abs([h.offset(1) ...
                         h.offset(1) + size(h.data,1)*h.binwidth(1)]));
    [ydata, kspike] = cross_intervals(t, ch, h.spikechan, maxinterv);
    if ~isempty(kspike),
        h.data = tdata_to_2Dhist(h, ydata, xdata(kspike(:,2)));
    end
    h.synccount = nsync;
    h.spikecount = count_events(ch,h.spikechan);
    
case 'PST-ShuffledCorrelation'
    [ydata, xdata] = shuffled_intervals(t, ch, h.syncchan, h.spikechan);
    if ~isempty(ydata),
        h.data = tdata_to_2Dhist(h, ydata, xdata);
    end
    h.synccount = count_events(ch,h.syncchan);
    h.spikecount = count_events(ch,h.spikechan);
    
case 'Period-Interval'
    [xdata, nperiod] = periodtimes(t, ch, h.period(index), h.syncchan);
    [ydata, kspike] = intervals(t, ch, h.spikechan, h.order);
    h.binwidth(2) = 1/size(h.data,2);
    h.data = tdata_to_2Dhist(h, ydata, xdata(kspike(:,2))/h.period(index));
    h.spikecount = length(kspike) + h.order;
    h.synccount = sum(ch == h.syncchan);
    h.periodcount = nperiod;
case 'Joint Interval'
    tdata = intervals(t, ch, h.spikechan, h.order);
    h.data = tdata_to_2Dhist(h, tdata(2:end), tdata(1:end-1));
    h.spikecount = length(tdata) + h.order;
otherwise
    error('Unsupported histogram type')
end


