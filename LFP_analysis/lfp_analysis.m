function analysis=lfp_analysis(lfp,frame_ts,deconv,fs,down_sample)
if nargin<5
    down_sample=400; % downsample to 300 Hz
end

factor=floor(fs/down_sample);
fs=fs/factor;

lfp=decimate(lfp,factor,'fir');

minute=floor(fs*60);
len=floor(length(lfp)/minute); %process by minute blocks
data=reshape(lfp(1:len*minute),minute,len);

wdw = 5*down_sample; %5s window
nol = floor(.5*wdw); %half overlapping windows
nfft = 2^16; 
f = 0:fs/nfft:fs/2;
%taking out 60Hz to normalize
fo1=55;
fi1=65; 
f_60out = f<fo1 | f>fi1;

nblock=size(data,2);
S=zeros(length(f),size(data,2));
for j=1:nblock
    data(:,j)=detrend(data(:,j)); %local detrend in minute block
    S(:,j) = pwelch(data(:,j),wdw,nol,nfft,fs);  
    S(:,j)=S(:,j)./sum(S(~f_60out,j)); %use 60hz as reference for power normalization
end
fo=.5;
fi=150;
fp = f>fo & f<fi;
t = 1:1:nblock;
figure
pcolor(t,f(fp),10*log10(S(fp,:)))
shading interp 
xlabel('time (min)')
ylabel('frequency (Hz)')

figure;
ax1=subplot(2, 1, 1);
fo=90;
fi=150;
fp = f>fo & f<fi;
t = 1:1:nblock;
pcolor(t,f(fp),10*log10(S(fp,:)))
shading interp 
ylabel('frequency (Hz)')

ax2=subplot(2, 1, 2);
fo=0.5;
fi=20;
fp = f>fo & f<fi;
t = 1:1:nblock;
pcolor(t,f(fp),10*log10(S(fp,:)))
shading interp 
linkaxes([ax1 ax2],'x');
xlabel('time (min)')
ylabel('frequency (Hz)')


%%%%%%%%%%
%detrend using 1s window with 500ms steps
lfp=locdetrend(lfp,fs,[1 .5]);
f60=designfilt('bandstopiir','FilterOrder',2,'HalfPowerFrequency1',55,'HalfPowerFrequency2',65,'DesignMethod','butter','SampleRate',fs);
lfp=filtfilt(f60,lfp);

f_delta=
f_theta=designfilt



% 
% % Define cycle
% RMS=squeeze(sqrt(mean(data.^2)));
% [H,N]=hist(RMS(nblock/2:end),100);   
% LIMinf=N(25); %low tresh 25%
% 
% LIMsup=N(75); %hi tresh 75%
% 
% figure;   
% 
% bar(N,H)
% hold on
% aux=1:max(H);
% plot(LIMinf,aux,'O-r');
% plot(LIMsup,aux,'o-r');
% 
% hold off;
% 
% figure
% % set(gcf,'color',[1 1 1],'Position',[  742   255   560   420]);  
% aux=ones(size(RMS(nblock/2:end),2));
% plot(RMS(nblock/2:end));
% hold on;
% plot(aux*LIMinf,'r--');
% plot(aux*LIMsup,'r');

    
% %% %%%%% ZCR %%%%%%
% 
% 
% for i=1:nblock
%     
% dado=datan(:,i);    
% ZCR=find(abs(diff(sign(dado)))==2);
% nZCR(i)=numel(ZCR);
% end
% 
% 
% 
% %%% SCATTER    
% figure
% hold on
% set(gcf,'color',[1 1 1],'Position',[ 35   255   560   420]); 
% set(gca,'FontName', 'Arial');
% scatter(nZCR(nblock/2:end),RMS(nblock/2:end),'k')
% box off;
% set(gca,'FontSize',18)
% ylabel('RMS')
% xlabel('nZCR')
% 
