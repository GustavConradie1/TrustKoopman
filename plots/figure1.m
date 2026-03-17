clear
%%
axissize=36;

%% plot pluto snapshots
load('states_pluto2.mat')
load('states_charon2.mat')
load('states_hydra2.mat')
load('states_styx2.mat')
load('states_nix2.mat')
load('states_kerberos2.mat')

DATA=[states_pluto2 states_charon2 states_hydra2 states_styx2 states_nix2 states_kerberos2];
DATA=DATA.';

%%
idx=1:1000;
figure
scatter3(DATA(7,idx)-DATA(1,idx),DATA(8,idx)-DATA(2,idx),DATA(9,idx)-DATA(3,idx),200,'.')
hold on
scatter3(DATA(13,idx)-DATA(1,idx),DATA(14,idx)-DATA(2,idx),DATA(15,idx)-DATA(3,idx),200,'.')
scatter3(DATA(25,idx)-DATA(1,idx),DATA(26,idx)-DATA(2,idx),DATA(27,idx)-DATA(3,idx),200,'.')
scatter3(DATA(31,idx)-DATA(1,idx),DATA(32,idx)-DATA(2,idx),DATA(33,idx)-DATA(3,idx),200,'.')
scatter3(DATA(19,idx)-DATA(1,idx),DATA(20,idx)-DATA(2,idx),DATA(21,idx)-DATA(3,idx),200,'.')
grid off
box off
axis off
view([-4.1857 21.8417])
[~,objh] = legend({'Charon','Hydra','Nix','Kerberos','Styx'},'interpreter','latex','fontsize',axissize,'location','best');
objhl = findobj(objh, 'type', 'patch');
set(objhl,'markersize', 40);
exportgraphics(gcf,'pluto_snapshots.png')

%%
clear
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
set(objhl, 'markersize', 10);
exportgraphics(gcf,'hyperparameter_optimization_pc.png','ContentType','vector','BackgroundColor','none')

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
exportgraphics(gcf,'hyperparameter_fixed_paths.png','ContentType','vector','BackgroundColor','none')

%% initial matern kernel
s=1/5; v=1;
range = linspace(10^(-16),100,10000);
vals = (s*range).^v.*besselk(v,s*range)/((s*range(1)).^v.*besselk(v,s*range(1)));

%% final matern kernel
s2 = 0.0522; v2=4.73;
vals2 = (s2*range).^v2.*besselk(v2,s2*range)/((s2*range(1)).^v2.*besselk(v2,s2*range(1)));

%% plot both
figure
plot(range,vals,'linewidth',linesize,'color',1.2*colors(1,:))
hold on
plot(range,vals2,'linewidth',linesize,'color',colors(2,:))
box on
ax=gca; ax.FontSize=axissize; axis([0 100 -0.001 1.001])
xlabel('$r$','interpreter','latex','fontsize',fontsize)
ylabel('$f(r)$','interpreter','latex','fontsize',fontsize)
legend({'Initial','Best'},'interpreter','latex','fontsize',fontsize,'location','best');
title('Radial basis function','interpreter','latex','fontsize',fontsize)
exportgraphics(gcf,'matern2.png')
