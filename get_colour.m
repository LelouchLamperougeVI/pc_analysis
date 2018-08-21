function c=get_colour(colour)
% don't like the built-in colormaps? Get a better one from here!

switch lower(colour)
    case 'black'
        c=[linspace(1,0,64)' linspace(1,0,64)' linspace(1,0,64)'];
    case 'red'
        c=[ones(64,1) linspace(1,0,64)' linspace(1,0,64)'];
    case 'blue'
        c=[linspace(1,0,64)' linspace(1,0,64)' ones(64,1)];
    case 'rainbow'
        c=rand(64,3);
    case 'magenta'
        c=[ones(64,1) linspace(1,0,64)' ones(64,1)];
end