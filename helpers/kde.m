function f = kde(x, xi, h, w, support)
% My own 2D kernel density estimation routine.
% Has two modes: First mode acts like a normal KDE, estimating the PDF at
% each xi query point. The second mode uses the kernel as weights to
% compute a weighted average for a given parameter. Second mode is
% activated by passing w as a vector (1 by default).
%
% If support is provided, apply boundary correction by reflection.
% support = [xlower ylower; xupper yupper];
%
% Examples
%
% temp = rand(1e4, 2);
% nbins = 100;
% xy = linspace(0, 1, nbins);
% [x, y] = meshgrid(xy, xy);
% f = kde(temp, [x(:), y(:)], .05);
% f = reshape(f, nbins, nbins);
% 
% imagesc(f)
% 
% f = kde(temp, [x(:), y(:)], .05, 1, [0, 0; 1, 1]);
% f = reshape(f, nbins, nbins);
% 
% imagesc(f)
% 
% f = kde(temp, [x(:), y(:)], .05, rand(length(temp), 1) .* 10, [0, 0; 1, 1]);
% f = reshape(f, nbins, nbins);
% 
% imagesc(f)

if ~ismatrix(x) || ~ismatrix(xi) || (size(x, 2) ~= 2) || (size(xi, 2) ~= 2)
    error('x and xi must be two columns matrices');
end

if nargin < 4
    w = 1;
end

if nargin > 4 % support provided, use reflection boundary correction
    support = mksupport(support);
else
    support = [nan nan];
end
support = permute(support, [3, 2, 1]);

f = nan(length(xi), 1);
norm = nan(length(xi), 1);

supported = support;
supported(isnan(supported)) = 0;
supported = supported - x;
supported = supported .* (-2 .* isnan(support) + 1);
for ii = 1:length(xi)
    d = supported - xi(ii, :);
    d = sum(d .^ 2, 2);
    d = d ./ (h^2) ./ 2;
    d = exp(-d);
    
    f(ii) = sum(d .* w, [1, 3]);
    norm(ii) = sum(d, [1, 3]);
end

if length(w) == 1
    norm = 2 * pi * h^2 * length(x);
end

f = f ./ norm;


function support = mksupport(support)
support = cat(1, [nan, nan], support);
support = support .* 2;
m = ...
    [1 1
     1 2
     1 3
     2 1
     2 2
     2 3
     3 1
     3 2
     3 3];
support = cat(2, support(m(:, 1), 1), support(m(:, 2), 2));