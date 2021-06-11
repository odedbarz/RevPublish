
clc
fignum = 11;
verbose = 1;

setup_environment('../');
   
drr   = get_DRR_list_and_indices;
n_drr = drr.n_drr;                  % # DRRs of used 
idx_dry = drr.ordered(1);




%% Load data
%   Run [main_aggregate_MUA_data.m] again to update this file if needed
% 
data_type   = 'SU';       % {'SU', MUA'}
fn.load.path= load.path_to_data('_data');
data_type   = upper(data_type);

switch data_type
    case 'SU'
        % Loads a struct with fields:
        %               H: [3600×6×150 double]
        %          S_list: {1×150 cell}
        %     neuron_list: [150×1 double]
        %         spec_st: [1×1 struct]
        %      tbl_impale: [437×20 table]
        fn.load.file = 'data_SU_(08-Jan-2021)_bw(5)_fbands(30)_win(NaN)ms_spec(gammatone).mat';       

    case 'MUA'
        %Loads a struct with fields:
        %               H: [7200×6×356 double]
        %        H_labels: [356×6 double]
        %     neuron_list: [356×1 double]
        %         spec_st: [1×1 struct]
        %         stim_st: [1×1 struct]
        %      tbl_impale: [437×20 table]        
        fn.load.file = 'data_MUA_(08-Jan-2021)_bw(5)_fbands(30)_win(NaN)ms_spec(gammatone).mat';        

    otherwise
        error('--> Unrecognized DATA_TYPE!');
        
end
fn.load.fullfile = fullfile( fn.load.path, fn.load.file );
data            = load(fn.load.fullfile);
n_units         = size(data.H,3); 
spec_st         = data.spec_st;
tbl_data        = data.(sprintf('tbl_%s', data_type));
duration_sec    = 36;      % (sec) 
assert(duration_sec == spec_st.duration_ms * 1e-3,...
    '--> ERROR: You are using the wrong stimulus duration!');

aux.vprint(verbose, '--> [main_loopover_stats.m] Loading file:\n\t...<%s>\n', fn.load.file);





%%
CC_units = nan(n_units, drr.n_drr);
H = data.H(:,drr.ordered,:);            % sort data by DRRs

for k = 1:n_units
    for rv = 2:drr.n_drr
        H1 = filter(hann(15),1,H(:,idx_dry,k));
        H2 = filter(hann(15),1,H(:,rv,k));
        CCk = corrcoef(H1, H2);

        CC_units(k,rv) = CCk(1,2);
    end
    
end


plot(CC_units(:,2), CC_units(:,end), '.');
axis square
xlim([0, 1]); ylim([0, 1]);
% histogram(CC_units(:,2), 25);
% histogram(CC_units(:,5), 25);
xlabel('CC ');
















