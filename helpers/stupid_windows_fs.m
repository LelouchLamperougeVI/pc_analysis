function deconv=stupid_windows_fs(deconv,frames_per_file)
% somewhere along the preprocessing, the stupid file system decided to
% order file numbers in the following way: 1, 10, 11, 12, ... 2, 3, 4, ...
% as a result, the deconv is all messed up
% pass deconv through this function to fix everything :)

% UPDATE: The culprit in the preprocessing has been identified and fixed.
% But this function is still relevant for older data

if nargin < 2
    frames_per_file=2000;
end

if size(deconv,1)/frames_per_file <= 9
    return
end

numFiles=size(deconv,1)/frames_per_file;

overflow = floor( (numFiles-9) / 10 ); % this is a hypothetical situation that hasn't occured yet, but just in case...
for ii = 1:overflow
    idx = (ii*frames_per_file + 1) : ((10+ii) * frames_per_file);
    buff = deconv( idx, : );
    deconv(idx, :) = [];
    deconv = [deconv; buff];
end

excess=mod( (numFiles-9), 10);
idx = ((overflow+1)*frames_per_file + 1) : ((excess + overflow + 1) * frames_per_file);

buff = deconv( idx, : );
deconv(idx, :) = [];
deconv = [deconv; buff];