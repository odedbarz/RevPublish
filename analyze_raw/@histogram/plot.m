function ph = plot(h, varargin)
% HISTOGRAM/PLOT - Plot histogram
% PLOT(H) plots the spike counts in histogram H (1-D or 2D).
% Plot uses bars for 1D histograms, and dot displays for 2D histograms.
%
% PLOT(H, 'rate') plots the histogram in units of discharge rate (spike/sec)
% rather than raw spike counts.
%
% For 1-D histograms, PLOT(H, 'line') plots H as lines rather than bars.
% For 2-D histograms, PLOT(H, 'image') plots H as an image.
% For dot rasters, PLOT(H, 'threshold', THRESHVAL) specifies a threshold value
% below which data are not plotted (default 0).
%
%  PLOT(H) returns a handle to line objects.
%
dmax = ndims(h.data);
if dmax > 3 || (dmax == 3 && length(h.value) > 1),
    plot4d(h, varargin{:});
    return;
    %    error('Can''t plot histograms/neurograms with more than 2 dimensions');
end

ph = [];    % Oded, 12/19/2018

how = 'bars';
what = 'count';
threshold = 0;
barWidth = 1;
if nargin < 1, varargin={}; end

while length(varargin) >= 1,
    switch lower(varargin{1})
        case 'rate'
            what = 'rate';
            h.data = h.data * rate_scale(h);
            varargin = varargin(2:end);
        case {'line', 'image', 'waterfall', 'bars', 'dots'}
            how = varargin{1};
            varargin = varargin(2:end);
        case 'threshold'
            if length(varargin) > 1,
                threshold = varargin{2};
                varargin = varargin(3:end);
            else
                varargin = varargin(2:end);
            end

        case 'barwidth',
            barWidth = varargin{2};
            varargin = varargin(3:end);

        otherwise
            break;
    end
end

if dmax == 3 && (strcmp(how, 'bars') || strcmp(how, 'line')),
    %     how = 'dots';
end

if ~isreal(h.data), h.data = abs(h.data); end

if dmax == 2 &&  min(size(h.data)) == 1,  % 1-D histogram

    dim = which_dim(h);
    xbins = get(h, 'BinTimes', dim);
    switch how
        case 'line'
            ph = plot(xbins, h.data, varargin{:});
        otherwise
            ph = bar(xbins + h.binwidth(dim)/2, h.data, barWidth, varargin{:});
            %         set(ph, 'FaceColor', 'b');
    end
    xaxis(xbins(1), xbins(end)+h.binwidth(dim));
    xlabel(axis_label(h, dim), 'FontWeight', 'bold');
    ylabel(axis_label(h, what), 'FontWeight', 'bold');


else               % 2-D histogram or neurogram

    if dmax == 2,   % 2-D histogram or neurogram of 1D histograms
        data = h.data;
        if length(h.value) == 1, % neurogram
            ybins = h.value{1};
        else                    % 2D histogram
            ybins = get(h, 'BinTimes', 1);
        end
    else            % neurogram of 2D histograms (dmax = 3)

%        xbins = get(h, 'BinTimes', dmax);
        yscale = BestScale(h.value{1});

        if isequal(how, 'dots'), % separate odd and even data to show in different colors
            %          dv = diff(h.value{1}(1:2));
            dv = 1;
            odd_data = [];
            yodd = [];
            even_data = [];
            yeven = [];
            for k = 1:size(h.data, 1),
                if mod(k, 2),
                    odd_data = [odd_data; squeeze(h.data(k,:,:))];
                    yodd = [yodd; k-1 + dv*[0:size(h.data,2)-1]'/size(h.data,2)];
                else
                    even_data = [even_data; squeeze(h.data(k,:,:))];
                    yeven = [yeven; k-1 + dv*[0:size(h.data,2)-1]'/size(h.data,2)];
                end
            end
            ybins(1) = min([yodd; yeven]);
            ybins(2) = max([yodd; yeven]);

        else
            % we collapse the second dimension (usually trial #) onto
            % the first (the neurogram variable).  First 2 dimensions must be
            % swapped so that Trial # varies fastest in collapsed dimension.
            data = permute(h.data, [2 1 3]);
            data = reshape(data, [size(data,1)*size(data,2) size(data,3)]);
            % interpolate neurogram variable to match length of collapsed dim
            ybins = interp1(0:size(h.data,1)-1, h.value{1}, ...
                [0:size(data,1)-1]/size(h.data,2), 'linear','extrap');
        end

    end
    xbins = get(h, 'BinTimes', dmax);

    switch how
        case 'image'   % plot as an image
            ph = imagesc(xbins, ybins, data);
            axis('xy');
        case 'waterfall'
            ph = waterfall(xbins, ybins', data);
            set(ph, 'facecolor', 'b');

        case {'line', 'bars'}, % added by keh to do PST neurograms as bar plots
            cla
            hold on;
            maxData = 1.05*max(max(data));
            if maxData <= 0, maxData = 1; end
            sync = get(h, 'SyncCount');
            for kp = 1:size(data,1),
                if sync(kp) > 0,
                    if isequal(how, 'line'),
                        ph = plot(xbins, data(kp, :)/maxData + kp-1, varargin{:});
                    elseif isequal(how, 'bars'),
                        ph = bar(xbins, data(kp, :)/maxData, barWidth, varargin{:});
                        set(ph, 'EdgeColor', 'b', 'FaceColor', 'b');
                        y = get(ph, 'YData') + kp-1;
                        set(ph, 'YData', y);
                        set(ph, 'BaseValue', kp-1);
                    end
                end
            end
            xaxis(xbins(1), xbins(end)+h.binwidth(2));

            yaxis(0, length(ybins)+1);
            yscale = BestScale(ybins);
            ytick = get(gca, 'YTick');
            if isequal(yscale, 'linear'),
                dy = min(diff(ybins));
                ytick = ytick*dy + ybins(1);
            else
                dy = min(log10(ybins(2:end)./ybins(1:end-1)));
                ytick = 10.^(ytick*dy)*ybins(1);
            end
            yticklabel = num2str(round(ytick(:))); % number of sig digits in labels needs to be thought about
            set(gca, 'YTickLabel', yticklabel);

            hold off;

        case 'dots',       % plot as dots
            cla;
            hold on;
            if dmax > 2,
                [ky, kx] = find(odd_data > threshold);
                [kye, kxe] = find(even_data > threshold);
                if ~isempty(yodd),
                    % Oded, 12/19/2018: ph --> ph(end+1)
                    ph(end+1) = plot(xbins(kx), yodd(ky), 'b.');
                end
                if ~isempty(yeven),
                    % Oded, 12/19/2018: ph --> ph(end+1)                    
                    ph(end+1) = plot(xbins(kxe), yeven(kye), 'm.');
                end
            else
                [ky, kx] = find(data > threshold);
                ph = plot(xbins(kx), ybins(ky), 'b.');
            end
            hold off;
            xaxis(xbins(1), xbins(end)+h.binwidth(2));
            yaxis(ybins(1), ybins(end));
            
    end    % switch how

    % MERGE - Works with neurograms but not single 2D histogram
    if ~isempty(h.value),
       set(gca, 'YTick', (1:size(h.data,1)) - 0.5);
       set(gca, 'YTickLabel', num2str(h.value{1}(:)));
    end
    xlabel(axis_label(h, dmax), 'FontWeight', 'bold');
    ylabel(axis_label(h, 1), 'FontWeight', 'bold');
    if strcmp(how, 'waterfall'),
        zlabel(axis_label(h, what), 'FontWeight', 'bold');
    end

end




