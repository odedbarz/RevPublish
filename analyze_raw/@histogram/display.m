function display(h)
% HISTOGRAM/DISPLAY - Display histogram parameters and data
%
disp(['Type:    ' h.type])
disp(['BinWidths: ' num2str(h.binwidth)])
disp(['Offsets: ' num2str(h.offset)])
if ~isempty(strfind(h.type, 'Period'))
    disp(['Period: ' num2str(h.period)])
end
disp(['SyncChan:   ' num2str(h.syncchan)])
disp(['SpikeChan:   ' num2str(h.spikechan)])

% neurogram variables
for k = 1:min(length(h.value),length(h.varname)),
    if length(h.value{k}) <= 12, 
         disp([h.varname{k} ':']);
         disp(h.value{k}');
    else disp(sprintf('%s: %d-element vector', ...
                      h.varname{k}, length(h.value{k})));
    end
end

if length(h.data) <= 8,
     disp(h.data)
else disp(sprintf('Data: %d X %d', size(h.data)))
end

disp(['SyncCount:   ' num2str(h.synccount(:)')])
if size(h.spikecount, 1) > 1 & size(h.spikecount, ndims(h.spikecount)) == 2,
    disp('SpikeCount:');
    disp(num2str(h.spikecount'));
else disp(['SpikeCount:  ', num2str(h.spikecount(:)')]);
end
    
