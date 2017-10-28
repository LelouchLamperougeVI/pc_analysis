% equally good at inducing an acid trip :)

freq=[2000 50];
x=1000;
y=100000;
nPlanes=200;
nChunks=200;

grat=@(theta,the,rho) rho.*cos(theta-the);
osc=@(p,t,theta) 0.5.*sin(2.*pi./p.*t+theta);

[the,rho]=cart2pol(repmat([1:x]',1,y),repmat(1:y,x,1));
% the=reshape(the,x,y);
% rho=reshape(rho,x,y);
% plane=zeros(x,y);
stack=[];
% the=gpuArray(the);
% rho=gpuArray(rho);

% plane=gpuArray(zeros(x,y));
% stack=gpuArray(zeros(x,y));
p=gpuArray(round(-diff(freq).*rand(1,nPlanes)+freq(2)));
theta=gpuArray(2*pi*rand(1,nPlanes));
ttheta=gpuArray(2*pi*rand(1,nPlanes));
rgb=rand(3);

for chunk=1:nChunks
    gPlane=gpuArray(zeros(x,y/nChunks));
    gStack=gpuArray(zeros(x,y/nChunks));
    gThe=gpuArray(the(:,(chunk-1)*y/nChunks+1:chunk*y/nChunks)); %(chunk-1)*y/nChunks*x+1:chunk*y/nChunks*x));
    gRho=gpuArray(rho(:,(chunk-1)*y/nChunks+1:chunk*y/nChunks)); %((chunk-1)*y/nChunks*x+1:chunk*y/nChunks*x));
    for i=1:nPlanes
        idx=grat(ttheta(i),gThe,gRho);
%         idx=reshape(idx,y/nChunks,x);
%         idx=idx';

        gPlane=osc(p(i),idx,theta(i));
        gStack=gStack+gPlane;
    end
%     gStack=(gStack-min(min(gStack)))./max(max(gStack));
%     gStack=gStack./max(max(gStack));
    if isempty(stack)
        stack=gather(gStack);
    else
        idx=gather(gStack);
        stack=[stack idx];
    end
end


stack=(stack-min(min(stack)))./max(max(stack));
stack=stack./max(max(stack));

stack=uint8(floor(stack.*255));
stack=ind2rgb(stack,prism(255));

imshow(stack);
% imagesc(stack)
% colormap gray
% pbaspect([y x 1]);
% set(gca, 'visible', 'off') ; 

figure;
imagesc(abs(log2(fftshift(fft2(rgb2gray(stack))))))
colormap gray
title('amplitude spectrum')
xlabel('cycles (pixels)')
ylabel('cycles (pixels)')
axis square

%%
imwrite(stack,'lol.tif')