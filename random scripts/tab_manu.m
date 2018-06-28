figure


subplot(1,3,1);
[f,x]=ecdf(tab_SI);
disp(['tab: ' num2str(mean(tab_SI)) ' +- ' num2str(sem(tab_SI))]);
plot(x,f,'r');
[f,x]=ecdf(tread_SI);
disp(['tread: ' num2str(mean(tread_SI)) ' +- ' num2str(sem(tread_SI))]);
hold on
plot(x,f,'k');
legend({'tablet','treadmill'},'location','southeast');
xlabel('mutual info. (bits)');
ylabel('cum. prob.');
[~,p]=ttest2(tab_SI,tread_SI)

subplot(1,3,2);
[f,x]=ecdf(tab_sparsity);
disp(['tab: ' num2str(mean(tab_sparsity)) ' +- ' num2str(sem(tab_sparsity))]);
plot(x,f,'r');
[f,x]=ecdf(tread_sparsity);
disp(['tread: ' num2str(mean(tread_sparsity)) ' +- ' num2str(sem(tread_sparsity))]);
hold on
plot(x,f,'k');
legend({'tablet','treadmill'},'location','southeast');
xlabel('sparsity (%)');
ylabel('cum. prob.');
[~,p]=ttest2(tab_sparsity,tread_sparsity)

subplot(1,3,3);
idx1=cellfun(@isempty, tab_width);
idx1=cell2mat(tab_width(~idx1)');
[f,x]=ecdf(idx1(:,1));
disp(['tab: ' num2str(mean(idx1(:,1))) ' +- ' num2str(sem(idx1(:,1)))]);
plot(x.*(100/80),f,'r');
idx2=cellfun(@isempty, tread_width);
idx2=cell2mat(tread_width(~idx2)');
[f,x]=ecdf(idx2(:,1));
disp(['tread: ' num2str(mean(idx2(:,1))) ' +- ' num2str(sem(idx2(:,1)))]);
hold on
plot(x.*(100/80),f,'k');
legend({'tablet','treadmill'},'location','southeast');
xlabel('place field width (cm)');
ylabel('cum. prob.');
[~,p]=ttest2(idx1(:,1),idx2(:,1))