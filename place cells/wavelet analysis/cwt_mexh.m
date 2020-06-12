function wt=cwt_mexh(signal,scale)
% Stupid matlab doesn't do cwt for mexican hat wavelet...
% Scale: upper limit for scale

if size(signal,1)~=1
    signal=signal';
end
if nargin<2
    scale=floor(length(signal)/2); % Nyquist freq thingie
end

wt=zeros(scale,length(signal));
for n=1:scale
    sd=n;
    kernel=ricker_wave(-.5*scale:.5*scale,sd);
    wt(n,:)=1/sqrt(sd).*conv(signal,kernel,'same');
end