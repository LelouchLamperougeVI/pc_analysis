function y=ricker_wave(x,sig)
% Ricker wavelet

y=2./(sqrt(3*sig)*pi^(1/4)).*(1-(x./sig).^2).*exp(-x.^2./2./sig^2);