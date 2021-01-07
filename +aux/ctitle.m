function title_str = ctitle(varargin)
%
%   function title_str = ctitle([ax], str1, str2);
%
% Centered title in Latex mode
%
% Using the template:
% '\begin{tabular}{c} The first line \\ Followed by another line \end{tabular}'
%

assert(2 <= length(varargin), '--> ERROR at [ctitle.m]: The input MUST contain at least TWO strings!!');

if isempty(varargin{1}) || ~ishandle(varargin{1}(1))
    ax = gca;
    argin_counter = 1;
else
    ax = varargin{1};
    assert(1 == length(ax), '--> [ERROR at ctitle]: AX must be ONE handle!!');
    argin_counter = 2;    
end

str1 = varargin{argin_counter};
str2 = varargin{argin_counter+1};

title_str = sprintf('\\begin{tabular}{c} %s \\\\ %s \\end{tabular}', str1, str2);

if 0 == nargout
    title(ax, sprintf('\\begin{tabular}{c} %s \\\\ %s \\end{tabular}', str1, str2));
end