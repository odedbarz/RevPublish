function t = histogram_types(which)

% HISTOGRAM_TYPES - List histogram types

% HISTOGRAM_TYPES with no argument returns the basic types

% HISTOGRAM_TYPES('derived') returns the derived types (e.g. synchrony)

% HISTOGRAM_TYPES('all') returns both basic and derived types.

%

t = {'none', 'PST', 'Interval', 'Period', 'Autocorrelation', 'Crosscorrelation', ...
     'Latency', 'Revcor', 'PST Raster', 'Interval Raster', 'Period Raster', ...
     'Autocorrelation Raster', 'Crosscorrelation Raster', 'PST-Interval', ...
     'PST-Period', 'PST-Autocorrelation', 'PST-Crosscorrelation', ...
     'Period-Interval', 'Joint Interval', 'PESE'};



if nargin > 0,

    dt = {'Synchrony Raster', 'PST-Synchrony', 'Synchrony-Interval'};



    switch which

    case 'derived'

        t = dt;

    case 'all'

        t = [t dt];

    end

end



    

    