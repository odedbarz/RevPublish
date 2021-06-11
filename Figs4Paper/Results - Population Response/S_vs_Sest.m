function [Tmu, Tse] = S_vs_Sest(CC, data_type, un, Tmu, Tse)
%
% function [Tmu, Tse] = S_vs_Sest(CC, data_type, un, Tmu, Tse)
%
% A simple Ad hoc fuction to extarct reconstruction data for all speakers from saved
% data table (Sdry vs. Sdrr or Sdrr vs. Sest).
%
% CC = data.tbl.CC.Variables;

% *** Sdry vs. Sdrr ***

% MAX_SU = 103;   % there are only 103 SUs in the data!
% if strcmpi('SU', data_type) && un > MAX_SU
%     col_name = sprintf('unit%d', un);
%     Tmu.([col_name, '_spk']) = nan(drr.n_drr, n_speakers);
%     return;     
% end

CCmu = squeeze( median(CC, 2) );
CCse = squeeze( mad(CC,[],2)/sqrt(size(CC,2)) );

col_name = sprintf('unit%d', un);
Tmu{:,col_name} = CCmu;
Tmu.([col_name, '_spk']) = CC;
Tmu.type = arrayfun(@(S) data_type, 1:size(CC,1), 'un', 0)';

Tse{:,col_name} = CCse;
Tse.type = arrayfun(@(S) data_type, 1:size(CC,1), 'un', 0)';
