%
% test_distcorr.m
%

clc;


data = load([load.path_to_data('data'), ...
    '/data_MUA_(08-Jan-2021)_bw(5)_fbands(30)_win(NaN)ms_spec(gammatone)_ALL-MEAS.mat'])



%%
drr = get_DRR_list_and_indices;
N = height(data.tbl_MUA_all);

dCorr = zeros(N, drr.n_drr);

for n = 1:N
    Tnm = data.tbl_MUA_all.all_meas{n};
    if 5 > height(Tnm)
        continue;
    end
    
    Tnm = Tnm(drr.ordered, :);
    dry_v = Tnm.single_meas{1};

    for m = 1:drr.n_drr
        drr_v = Tnm.single_meas{m};
        n_meas = min( [size(dry_v,2), size(drr_v,2)] );
        
        dCorr(n,m) = distcorr.distcorr( dry_v(:,1:n_meas), drr_v(:,1:n_meas) );
        
    end
   
    disp( dCorr(n,:) )
    
end






