function [Y, fs] = load_wav_files(file_list)
%
% function Y = load_files(files)
%
%

len_files = length(file_list);
Y = cell(1,len_files);

fs = [];

for jj = 1:len_files
    [Y{jj}, fs_jj] = audioread( file_list{jj} );
    if isempty(fs)
        fs = fs_jj;
    elseif fs ~= fs_jj
        error('--> ERROR at [load_wav_files.m]: you are loading WAV files with different sampling frequencies!!');
    end
    
    
end







