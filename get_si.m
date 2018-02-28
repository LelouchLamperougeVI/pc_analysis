function SI=get_si(raw_psth)

temp=[sum(raw_psth>0);sum(raw_psth==0)];
temp=temp./sum(sum(temp));
temp=temp.*log2(temp./sum(temp,1)./sum(temp,2));
temp(isnan(temp))=0;
SI=sum(sum(temp));
SI=reshape(SI,1,[]);