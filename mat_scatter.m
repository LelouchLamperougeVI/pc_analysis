function h=mat_scatter(data,m_size)
% Create scatter plot from matrix

if nargin<2
    m_size=1;
end

[m,n]=find(data);
h=plot(m,n,'.','markersize',m_size);