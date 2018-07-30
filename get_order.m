function order=get_order(analysis)

order=analysis.stack;
[~,order]=max(order);
[~,order]=sort(order);