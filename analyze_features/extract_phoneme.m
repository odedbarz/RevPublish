function [tbl_phn, tbl_metadata] = extract_phoneme(phn, binwidth, duration_sec)
%
%   function [tbl_phn, tbl_metadata] = extract_phoneme(phn, binwidth, duration_sec)
%
% Inputs:
%   phn
%   duration_sec

if nargin < 3
    duration_sec = 36;  % (sec)
end


%
fn.tbl.path     = load.path_to_data('Stimulus');
fn.tbl.file     = sprintf('metadata_(%d)_wav_(30-Jun-2020)', duration_sec);
dummy           = load( fullfile( fn.tbl.path, fn.tbl.file ) );
tbl_metadata    = dummy.tbl_metadata;

fs = 1/(1e-3*binwidth);     % (Hz)

% phn = 's';
var_names = {'phn', 'n0', 'n1', 't0', 't1'};
tbl_phn = table( {''}, 'VariableNames', {'phn'});
tbl_phn = [tbl_phn, array2table( nan(1,length(var_names)-1), 'VariableNames', var_names(2:end))];

row = 0;    % last table row
for k = 1:height(tbl_metadata)
    phn_loc = find( strcmpi(tbl_metadata.phn{k}, phn) );

    % starting time of the k'th speaker
    n0_sp = tbl_metadata.t0(k);
    fs_wav = tbl_metadata.fs(k);
    t0_sp = n0_sp/fs_wav;
        
    % Add the phoneme to the table
    for jj = 1:numel(phn_loc)
        phn_loc_j = phn_loc(jj);
                
        % make sure that it works properly
        assert( isequal( tbl_metadata.phn{k}{phn_loc_j}, phn), '--> Somthing is WRONG!' );   
        row = row + 1;
        warning off
        tbl_phn.phn{row} = tbl_metadata.phn{k}{phn_loc_j};
        warning on
        
        % start & end times, in seconds
        tbl_phn.t0(row) = t0_sp + 1/fs_wav * tbl_metadata.phn_start{k}(phn_loc_j);
        tbl_phn.t1(row) = t0_sp + 1/fs_wav * tbl_metadata.phn_end{k}(phn_loc_j);
        
        % start & end times, in samples
        tbl_phn.n0(row) = floor(tbl_phn.t0(row)*fs);
        tbl_phn.n1(row) = floor(tbl_phn.t1(row)*fs);
    end
    
end


assert( 0 < row, sprintf('ERROR: could not find the phoneme %s in the database', phn) );


