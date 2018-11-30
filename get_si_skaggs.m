function SI=get_si_skaggs(fr,mu_fr,Pi)

Pi=Pi./sum(Pi);

SI=Pi.*fr./mu_fr.*log2(fr./mu_fr);
SI(isnan(SI))=0;
SI=sum(SI);