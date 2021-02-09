function CCs = CC_stimuli( spec_st )
%
% function CCs = CC_stimuli( spec_st )
%
% Description:
% Computes correlation coefficients (CCs) between STIMULI.
%
%

drr = get_DRR_list_and_indices;

Sdry = spec_st.Sft{drr.dry};
CCs = nan(1, drr.n_drr);    % CCs of responses
for k = 1:drr.n_drr
    rvi = drr.ordered(k);
    Sk = spec_st.Sft{rvi};
    
    % CCs: correlation between SRY & DRR stimuli    
    gof = goodness(Sdry, Sk);
    CCs(k) = gof.CC;
end

