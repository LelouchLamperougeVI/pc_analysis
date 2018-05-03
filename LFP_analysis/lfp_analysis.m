function analysis=lfp_analysis(lfp,frame_ts,deconv,fs,down_sample)
if nargin<5
    down_sample=300; % downsample to 300 Hz
end

factor=floor(fs/down_sample);
fs=fs/factor;

lfp=decimate(lfp,factor,'fir');

minute=floor(fs*60);
len=floor(length(lfp)/minute); %process by minute blocks
data=reshape(lfp(1:len*minute),minute,len);

nblock=size(data,2);
for j=1:nblock
    wdw = 5*fs; nol = 0.5*wdw; nfft = 2^16; 
    f = [0:fs/nfft:fs/2];
    %taking out 60Hz to normalize
    fo1=55;
    fi1=65; 
    f_60out = f<fo1 | f<fi1;
    datan(:,j)=detrend(data(:,j));
    [S(:,j),w] = pwelch(data(:,j),wdw,nol,nfft,fs);  
    Sn(:,j)=S(:,j)./sum(S(f_60out,j));
end
fo=.5;
fi=100;
fp = f>fo & f<fi;
t = 1:1:nblock;
figure
pcolor(t,f(fp),10*log10(Sn(fp,:)))
shading interp 

figure;
ax1=subplot(2, 1, 1);
fo=100;
fi=200;
fp = f>fo & f<fi;
t = 1:1:nblock;
pcolor(t,f(fp),10*log10(Sn(fp,:)))
shading interp 

ax2=subplot(2, 1, 2);
fo=0.5;
fi=20;
fp = f>fo & f<fi;
t = 1:1:nblock;
pcolor(t,f(fp),10*log10(Sn(fp,:)))
shading interp 
linkaxes([ax1 ax2],'x');



%% %define cicle

RMS=squeeze(sqrt(mean(datan.^2,1)));
[H,N]=hist(RMS(nblock/2:end),100);   
LIMinf=N(25); %low tresh 25%

LIMsup=N(75); %hi tresh 75%

figure;   

bar(N,H)
hold on
aux=1:max(H);
plot(LIMinf,aux,'O-r');
plot(LIMsup,aux,'o-r');

hold off;

figure
% set(gcf,'color',[1 1 1],'Position',[  742   255   560   420]);  
    aux=ones(size(RMS(nblock/2:end),2));
    plot(RMS(nblock/2:end));
    hold on;
    plot(aux*LIMinf,'r--');
    plot(aux*LIMsup,'r');
        
    
%% %%%%% ZCR %%%%%%


for i=1:nblock
    
dado=datan(:,i);    
ZCR=find(abs(diff(sign(dado)))==2);
nZCR(i)=numel(ZCR);
end



%%% SCATTER    
figure
hold on
set(gcf,'color',[1 1 1],'Position',[ 35   255   560   420]); 
set(gca,'FontName', 'Arial');
scatter(nZCR(nblock/2:end),RMS(nblock/2:end),'k')
box off;
set(gca,'FontSize',18)
ylabel('RMS')
xlabel('nZCR')

