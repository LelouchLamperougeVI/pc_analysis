function dirTree(s, offset)
% recursively print the tree structure of directory
% s             struct fith fields 'name' and 'child'
%   s.name      char vector
%   s.child     same as s (recursive)

if nargin < 2
    offset = 0;
end
if isempty(s)
    return
end

clen = {s.name};
clen = ceil( max( cellfun(@length, clen) ) / 2 );

for ii = 1:length(s)
    disp(['|' repmat(' ', 1, offset) '|'])
    disp(['|' repmat(' ', 1, offset) '|'])
    disp(['|' repmat(' ', 1, offset) '|--> ' s(ii).name ])
    
    if isfield(s, 'child')
        dirTree( s(ii).child, offset + clen + 4 );
    end
end
