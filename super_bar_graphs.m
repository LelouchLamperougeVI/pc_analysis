function h = super_bar_graphs(means,error,groups,pval)

h=bar(means);
hold on;
errorbar(means,error,'k.','linewidth',1);
sigstar(groups,pval);