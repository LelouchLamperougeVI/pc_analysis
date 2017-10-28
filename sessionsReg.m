function dd = sessionsReg(moving,fixed)
% Register mean images across recordings
%   Inputs:
%       moving & fixed: separate mean images
%   Outputs:
%       dd: delta for transform

fixed=imread(fixed);
moving=imread(moving);
subplot(2,1,1);
imshowpair(fixed,moving,'montage')
dd=imregcorr(moving,fixed);
Rfixed = imref2d(size(fixed));
movingReg = imwarp(moving,dd,'OutputView',Rfixed);
subplot(2,1,2);
imshowpair(fixed,movingReg,'montage');