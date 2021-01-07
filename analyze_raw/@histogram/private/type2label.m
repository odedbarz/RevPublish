 function s = type2label(type)
 % TYPE2LABEL - Convert histogram type to axis label
 %
 switch type
    case 'PST'
        s = 'Peristimulus Time (ms)';
    case 'Interval'
        s = 'Interspike Interval (ms)';
    case 'Period'
        s = 'Period Phase (cycles)';
    case 'Autocorrelation' 
        s = 'All-Order Interval (ms)';
    case 'Crosscorrelation' 
        s = 'Cross Interval (ms)';
    case 'Recurrence'
        s = 'Recurrence Time (ms)';
    case 'ShuffledCorrelation'
        s = 'Shuffled Interval (ms)';
    case 'Joint'
        s = 'Preceding Interval (ms)';
    case 'Latency'
        s = 'Latency (ms)';
    case {'Revcor', 'PESE'}
        s = 'Pre-Spike Time (ms)';
    case 'Raster'
        s = 'Stimulus Trial';
    otherwise
        s = 'Time (ms)';
end
