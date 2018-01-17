function h=plot_sequences(sequences)
sequences(sequences==0)=NaN;
[~,order]=sort(mean(sequences,2));
sequences=sequences(order,:);
m=max(sequences)+6;
m=cumsum(m)-m(1);
linear=sequences+m;
linear=reshape(linear,1,[]);
y=repmat(1:size(sequences,1),1,size(sequences,2));
h=plot(linear,y,'.');
hold on
stem(m(2:end)-3,ones(1,size(sequences,2)-1)*size(sequences,1),'marker','none');