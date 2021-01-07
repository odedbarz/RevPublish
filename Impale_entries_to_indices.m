function midx = Impale_entries_to_indices(S)
%
%   function midx = Impale_entries_to_indices(S)
%
% Input:
% S: (1x1) Impale structure or DRR table. 
%
% Description: 
% Use this function to extract the available DRR measurements and set them
% into a 1x6 long index vector.
%
%     #         1        2       3        4        5       6
%     DRR      Dry  9.4 dB  4.8 dB  -2.5 dB  -8.2 dB      --
%     Revb.   100%     80%     80%      20%      20%    100%           
%     Dist.   1.5m    1.5m    3.0m     1.5m     3.0m    3.0m

drr = get_DRR_list_and_indices;
n_drr = length(drr.labels);
midx = false(1, n_drr);

if isstruct(S)
    %{
    % Impale structure
    values_dist = S.outerSeq.master.values;    % Source distance
    
    % Legacy, set all the distances & reverberation coefficients to have
    % the same format
    if 0.1*values_dist(1) > 10
        values_dist   = 0.01*values_dist;  % (m) distance between source and "ears"
    else
        values_dist   = 0.1*values_dist;   % (m) distance between source and "ears"
    end
    
    values_reverb = S.innerSeq.master.values;    % Walls absorption coefficients    
    if 10*values_reverb(1) <= 100
        %values_reverb = values_reverb;       % (percent) absorption coefficient
    %else
        values_reverb = 10*values_reverb;    % (percent) absorption coefficient
    end

    for ii = 1:length(values_dist)
        dist_ii = values_dist(ii);

        for jj = 1:length(values_reverb)
            reverb_jj = values_reverb(jj);      % (percent) absorption coefficient
            idx_ij = (drr.dist == dist_ii) & (drr.revb == reverb_jj);
            assert(~isempty(idx_ij), '--> [Impale_entries_to_indices.m]: something is WRONG with this measurement!');

            %c = c+1;
            midx(idx_ij) = true;

        end
    end   
     %}
    
    nsync = cellfun(@(CH) nnz(0==CH), S.ch, 'UniformOutput', true);            
    midx(1:length(nsync(:))) = nsync(:) > 0;
    
elseif istable(S)
    % MATLAB's table (as in the DRR table)
    values_dist = S.Dist;           % Source distance
    values_reverb = S.Reverb;       % Walls absorption coefficients
    
    for ii = 1:length(values_dist)
        dist_ii = values_dist(ii);
        reverb_ii = values_reverb(ii);      % (percent) absorption coefficient
        idx_ii = (drr.dist == dist_ii) & (drr.revb == reverb_ii);

        if ~any(idx_ii)
            continue;
        end

        %c = c+1;
        midx(idx_ii) = true;        
    end
    
else
    error('--> [Impale_entries_to_indices.m]: unrecognized input!');
end


% drr = get_DRR_list_and_indices;
% n_drr = length(drr.labels);
% midx = false(1, n_drr);
% %c = 0;

% for ii = 1:length(values_dist)
%     dist_ii = values_dist(ii);
%     
%     for jj = 1:length(values_reverb)
%         reverb_jj = values_reverb(jj);      % (percent) absorption coefficient
%         idx_ij = (drr.dist == dist_ii) & (drr.revb == reverb_jj);
%         assert(~isempty(idx_ij), '--> [Impale_entries_to_indices.m]: something is WRONG with this measurement!');
%         
%         %c = c+1;
%         midx(idx_ij) = true;
%         
%     end
% end

% % for ii = 1:length(values_dist)
% %     dist_ii = values_dist(ii);
% %     reverb_ii = values_reverb(ii);      % (percent) absorption coefficient
% %     idx_ii = (drr.dist == dist_ii) & (drr.revb == reverb_ii);
% % 
% %     if ~any(idx_ii)
% %         continue;
% %     end
% % 
% %     %c = c+1;
% %     midx(idx_ii) = true;        
% % 
% % end


