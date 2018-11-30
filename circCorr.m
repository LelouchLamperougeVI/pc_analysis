function h=circCorr(r,D,assemblies,thres,no_order)
% draw a circular connectivity graph between ensembles

if nargin<5
    no_order=false;
end

n=length(r);
angles=0:2*pi/n:2*pi;
[x,y]=pol2cart(angles,ones(1,n+1));
x(end)=[];y(end)=[];
angles(end)=[];
angles=rad2deg(angles);

if no_order
    all_neurons=1:n;
else
    all_neurons=cell2mat(assemblies);
    all_neurons=[all_neurons setxor(all_neurons,1:n)];
end

h=figure;
hold on
set(gca,'Visible','off');
axis square

count=0;
for i=1:length(assemblies)
    if no_order
        l=plot(x(assemblies{i}),y(assemblies{i}),'o');
    else
        l=plot(x(1+count:length(assemblies{i})+count),y(1+count:length(assemblies{i})+count),'o');
    end
    color=l.Color;
    for j=1:length(assemblies{i})
        count=count+1;
        if no_order
            text(x(assemblies{i}(j)),y(assemblies{i}(j)),num2str(all_neurons(assemblies{i}(j))),'Rotation',angles(assemblies{i}(j)),'Color',color);
        else
            text(x(count),y(count),num2str(all_neurons(count)),'Rotation',angles(count),'Color',color);
        end
    end
end

if no_order
    idx=setxor(cell2mat(assemblies),1:n);
    plot(x(idx),y(idx),'ko');
    for j=1:length(idx)
        count=count+1;
        text(x(idx(j)),y(idx(j)),num2str(idx(j)),'Rotation',angles(idx(j)),'Color',zeros(1,3));
    end
else
    plot(x(1+count:end),y(1+count:end),'ko');
    for j=1:n-length(cell2mat(assemblies))
        count=count+1;
        text(x(count),y(count),num2str(all_neurons(count)),'Rotation',angles(count),'Color',zeros(1,3));
    end
end

mask=true(n);
mask(cell2mat(assemblies),cell2mat(assemblies))=false;
mask=mask | D>thres;
r=r.*~diag(ones(1,n));
r=squareform(r);
color=linspace(min(r),max(r),64);
opac=linspace(0,1,64);
r=knnsearch(color',r');
r=squareform(r);
r(mask)=nan;

color=hot;
for i=1:n
    for j=i+1:n
        if ~isnan(r(i,j))
            if no_order
                l=plot(x([i j]),y([i j]),'color',color(r(i,j),:));
            else
                l=plot(x([find(all_neurons==i) find(all_neurons==j)]),y([find(all_neurons==i) find(all_neurons==j)]),'color',color(r(i,j),:));
            end
            
            alpha(l,opac(r(i,j)));
        end
    end
end
