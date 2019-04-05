function I = generate_gratings(varargin)
% generate grayscale gratings
% Inputs:
%   'dims':         dimensions of image (default - [96 96])
%   'f':            frequency of oscillation (default - .1 cycles/pixel)
%   'phase':        phase in degrees (default - 0)
%   'orientation':  orientation in degrees (default - 0)

[dims, f, psi, orientation] = parse_inputs(varargin);

[x,y] = ind2sub(dims, 1:prod(dims));
[theta,rho] = cart2pol(x,y);

proj = abs(rho).*cos(theta-deg2rad(orientation));
I = .5 .* sin(2*pi * f * proj + psi);
I = reshape(I, dims) + .5;


function [dims, f, psi, orientation] = parse_inputs(inputs)

dims = [96 96];
f = .1;
psi = 0;
orientation = 0;

count = 1;
while count < length(inputs)
    switch lower(inputs{count})
        
        case 'dims'
            dims = inputs{count+1};
        
        case 'f'
            f = inputs{count+1};
            
        case 'orientation'
            orientation = inputs{count+1};
            
        case 'phase'
            psi = deg2rad(inputs{count+1});
            
    end
    count = count+2;
end