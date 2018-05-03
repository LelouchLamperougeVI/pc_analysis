function d=bi2de(b)
% fuck you matlab, I'm not installing your stupid toolbox

% if size(b,1)==1
%     b=b';
% end

d=(2.*ones(1,size(b,1))).^(0:size(b,1)-1);

if size(b,1)~=1
    d=sum(b.*d(end:-1:1)');
else
    d=b.*d;
end