function mean2tiff(in,out)
% Turn mean images to tiff
load(in);
mimg=(mimg-min(min(mimg)))./max(max(mimg));
fig=figure;
h=imshow(mimg);
imwrite(double(h.CData),out);
close(fig);