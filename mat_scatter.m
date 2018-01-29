function h=mat_scatter(data,param)
% Create scatter plot from matrix

ops='.';

if nargin>1
    ops=[ops param];
end

[m,n]=find(data);
h=plot(m,n,ops,'markersize',1);