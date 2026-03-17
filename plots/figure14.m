clear
load('optimization_results_good.mat')
load('hyperparameter_fixed_results.mat')
fontsize=36;
axissize=28;
linesize=5;
colors=orderedcolors('gem');
purple = [150/256 50/256 150/256];
worst_loss=0.0858;

%% hyperparameter optimization over epochs
green=2.2*[11/256 64/256 0/256];
figure
semilogy(0:1:500,exact_error_history,'linewidth',linesize,'color','black')
hold on
plot(0:1:500,loss_history,'linewidth',linesize,'color',colors(1,:))  
plot(0:1:500,pointwise_error_history,'linewidth',linesize,'color',green)  
scatter(1,worst_loss,1500,purple,'x','linewidth',linesize)
scatter(best_epoch,best_loss,1500,purple,'*','linewidth',linesize)
text(45,worst_loss-0.01,'Initial epoch','fontsize',28,'interpreter','latex','color',purple)
text(40,worst_loss-0.025,'$s=0.2,\nu=1$','fontsize',28,'interpreter','latex','color',purple)
text(best_epoch-55,best_loss-0.003,'Best epoch','fontsize',28,'interpreter','latex','color',purple)
text(best_epoch-95,best_loss-0.006,'$s=0.0522,\nu=4.73$','fontsize',28,'interpreter','latex','color',purple)
ax=gca; ax.FontSize=axissize;
title('Pluto-Charon system','fontsize',fontsize,'interpreter','latex')
xlabel('Epoch','interpreter','latex','fontsize',fontsize)
ylabel('Error','interpreter','latex','fontsize',fontsize)
xticks([0,100,200,300,400,500])
yticks([10^(-2) 10^(-1)])
box on
[~,objh] = legend({'Exact $L^2$ error','$L^2$ error bound','Exact pointwise error$\;\;$'},'interpreter','latex','fontsize',24,'location','best');
objhl = findobj(objh, 'type', 'patch');
set(objhl,'markersize', 10);
exportgraphics(gcf,'hyperparameter_optimization_pc.pdf','ContentType','vector','BackgroundColor','none')

%% comparison of errors for best and worst results 
figure
semilogy(0:1:20,exact_errors_init,'linewidth',linesize,'color','black')
hold on
plot(0:1:20,error_bounds_init,'linewidth',linesize,'color',1.2*colors(1,:))  
plot(0:1:20,exact_errors_best,':','linewidth',linesize,'color','black')  
plot(0:1:20,error_bounds_best,':','linewidth',linesize,'color',1.2*colors(1,:))  
scatter(0:5:20,error_bounds_init(1:5:end),500,purple,'x','linewidth',linesize)
scatter(0:5:20,error_bounds_best(1:5:end),800,purple,'*','linewidth',linesize)
box on
ax=gca; ax.FontSize=axissize; axis tight
title('Pluto-Charon system','fontsize',fontsize,'interpreter','latex')
xlabel('Time (days)','interpreter','latex','fontsize',fontsize)
ylabel('$L^2$ Error','interpreter','latex','fontsize',fontsize)
legend({'Exact, initial$\;\;$','Bound, initial$\;\;$','Exact, best','Bound, best'},'interpreter','latex','fontsize',24,'location','southeast')
exportgraphics(gcf,'hyperparameter_fixed_paths.pdf','ContentType','vector','BackgroundColor','none')
