function ph = plot4d(h, varargin)
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
% if dmax ~= 4 || length(h.value) ~= 2,
%    error('Can''t plot histograms/neurograms with more than 2 dimensions');
% end

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

if ~isreal(h.data), h.data = abs(h.data); end

if isequal(how, 'dots'),
   dy = (0:size(h.data, 3) - 1)' / size(h.data, 3);
   dx = (0:size(h.data, 4)- 1)' / size(h.data, 4);
   
   odd_data = [];
   xodd = [];
   yodd = [];
   even_data = [];
   xeven = [];
   yeven = [];
   
   for k = 1:size(h.data, 1),
      for j = 1:size(h.data, 2),
         if mod(k+j, 2),
            odd_data = [odd_data; squeeze(h.data(k,j,:,:))];
            xodd = [xodd; (k-1)*ones(size(dy))];
            yodd = [yodd; (j-1) + dy(:)];
         else
            even_data = [even_data; squeeze(h.data(k,j,:,:))];
            xeven = [xeven; (k-1)*ones(size(dy))];
            yeven = [yeven; (j-1) + dy(:)];
         end
      end
   end
   xbins(1) = min([xodd; xeven]);
   xbins(2) = max([xodd; xeven]);
   ybins(1) = min([yodd; yeven]);
   ybins(2) = max([yodd; yeven]);

else
   % we collapse the second dimension (usually trial #) onto
   % the first (the neurogram variable).  First 2 dimensions must be
   % swapped so that Trial # varies fastest in collapsed dimension.
   dx = (0:size(h.data, 3)- 1)' / size(h.data, 3);
   data = h.data;
   xbins = [0 size(h.data, 1)];
   ybins = [0 size(h.data, 2)];
%    data = permute(h.data, [2 1 3]);
%    data = reshape(data, [size(data,1)*size(data,2) size(data,3)]);
%    % interpolate neurogram variable to match length of collapsed dim
%    ybins = interp1(0:size(h.data,1)-1, h.value{1}, ...
%       [0:size(data,1)-1]/size(h.data,2), 'linear','extrap');
end

% xbins = get(h, 'BinTimes', dmax);

switch how
%    case 'image'   % plot as an image
%       ph = imagesc(xbins, ybins, data);
%       axis('xy');
%    case 'waterfall'
%       ph = waterfall(xbins, ybins', data);
%       set(ph, 'facecolor', 'b');

   case {'line', 'bars'}, % added by keh to do PST neurograms as bar plots
      hold on;
      maxData = 1.05*max(data(:));
      if maxData <= 0, maxData = 1; end
      sync = get(h, 'SyncCount');
      for kp = 1:size(data, 1),
         for jp = 1:size(data, 2),
            if sync(kp, jp) > 0,
               if isequal(how, 'line'),
                  ph = plot(xbins, data(kp, :)/maxData + kp-1, varargin{:});
               elseif isequal(how, 'bars'),
                  ph = bar((kp-1)+dx, squeeze(data(kp, jp, :))/maxData, barWidth, varargin{:});
                  set(ph, 'EdgeColor', 'b', 'FaceColor', 'b');
                  y = get(ph, 'YData') + jp - 1;
                  x = get(ph, 'XData');
                  set(ph, 'XData', x, 'YData', y, 'BaseValue', jp-1);
               end
            end
         end
      end
      xaxis(xbins(1), xbins(end));
      yaxis(ybins(1), ybins(end));

   case 'dots',       % plot as dots
      cla;
      set(gca, 'TickDir', 'out');
      hold on;
      [ky, kx] = find(odd_data > threshold);
      [kye, kxe] = find(even_data > threshold);
      if ~isempty(yodd),
         ph = plot(xodd(ky)+dx(kx), yodd(ky), 'b.');
      end
      if ~isempty(yeven),
         ph = plot(xeven(kye)+dx(kxe), yeven(kye), 'm.');
      end
%       hold off;
%       xaxis(xbins(1), xbins(end)+h.binwidth(2));
      xaxis(xbins(1), xbins(end)+1);
      yaxis(ybins(1), ybins(end));

end

xlabel(axis_label(h, 1), 'FontWeight', 'bold');
ylabel(axis_label(h, 2), 'FontWeight', 'bold');
set(gca, 'XTick', 0.5 + (0:size(h.data, 1)));
set(gca, 'XTickLabel', num2str(h.value{1}));
set(gca, 'YTick', 0.5 + (0:size(h.data, 2)));
set(gca, 'YTickLabel', num2str(h.value{2}));

if isequal(how, 'dots'),
   yl = get(gca, 'YLim');
   for x = 1:size(h.data,1)-1,
      plot(x*[1 1],yl,'k-');
   end
   hl = get(gca, 'XLabel');
   xlabPos = get(hl, 'Position');

%    [ax1,ay1]=dsxy2figxy(gca,0,xlabPos(2));
%    [ax2,ay2]=dsxy2figxy(gca,1,xlabPos(2));

%    tmax = size(h.data,4)*h.binwidth(2);
%    ha = annotation('textarrow', [ax1 ax2], [ay1 ay2]);
%    set(ha, 'String', sprintf('%d ms',round(tmax)));
end





