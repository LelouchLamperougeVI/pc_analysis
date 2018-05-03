function [sce,H]=entropy_sce(deconv,scale,plotFlag)
% Gradient ascent over "entropy surface" to detect SCEs
% Set 2^scale ~= maximum number of frames for (re-)activation window

if nargin<3
    plotFlag=0;
end

if islogical(deconv)
    deconv=double(deconv);
end

H=zeros(size(deconv,1),scale+1);
for i=0:scale
%     H(:,i+1)=loc_entropy(deconv,2^i,5);
    H(:,i+1)=bin_loc_entropy(deconv,2^i);
end
H(isnan(H))=0;
%null dist
H_null=zeros(size(deconv,1),scale+1);
shuffled_d=mat_circshift(deconv,randi(size(deconv,1),1,size(deconv,2)));
for i=0:scale
    H_null(:,i+1)=bin_loc_entropy(shuffled_d,2^i);
end
H_null(isnan(H_null))=0;
H_null=mean(H_null);

H=H./H_null;

H=fast_smooth(H,1.5);

sce=false(size(H));
sce(:,end)=imregionalmax(H(:,end));
% sce(:,end)=H(:,end);
for i=size(H,2)-1:-1:1
%     maxi=imregionalmax(H(:,i));
%     maxi=H(:,i).*maxi;
    maxi=H(:,i);
%     region_maxi=H(:,i);
    for j=find(sce(:,i+1))'
        temp=maxi;
%         temp2=region_maxi;
        try
            temp(1:j-floor(2^i/4)-1)=0;
%             temp2(1:j-floor(2^i/4)-1)=0;
        catch
            disp('jinx')
        end
        try
            temp(j+floor(2^i/4)+1:end)=0;
%             temp2(j+floor(2^i/4)+1:end)=0;
        catch
            disp('jinx')
        end
        temp=imregionalmax(temp).*temp;
        temp=temp>H(j,i+1);
        sce(:,i)=sce(:,i)|temp;
        if any(temp)
            sce(j,i+1)=false;
        end
    end
end

if plotFlag
    figure;
    ax1=subplot(2,1,1);
    imagesc(fast_smooth(deconv,5)');
    ax2=subplot(2,1,2);
    imagesc(H');
    hold on
    [x,y]=find(sce);
    plot(x,y,'*r');
    linkaxes([ax1 ax2],'x');
end
