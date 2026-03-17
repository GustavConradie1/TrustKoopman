%%
fontsize=36;
axissize=28;
linesize=5;
colors=orderedcolors('gem');

%%
load('duffing_kernel_error_bounds.mat')
normalization=sqrt(normalization);

%% Plot relative forecast errors compared to test data
figure
semilogy(0:1:steps,exact_errors/normalization,'linewidth',linesize,'color','black')
hold on
plot(0:1:steps,(min(proj_errors,proj_errors2)+delta_errors)/normalization,'linewidth',linesize,'color',1.2*colors(1,:))    
box on
ax=gca; ax.FontSize=axissize;
title('Duffing oscillator','fontsize',fontsize,'interpreter','latex')
xlabel('Time steps','interpreter','latex','fontsize',fontsize)
ylabel('RKHS Error','interpreter','latex','fontsize',fontsize)
legend({'Exact','Bound'},'interpreter','latex','fontsize',fontsize,'location','northwest')
exportgraphics(gcf,'duffing_oscillator_kernel_error_bound.pdf','ContentType','vector','BackgroundColor','none')

%%
load('lorenz_kernel_error_bounds.mat')

%% 
figure
semilogy((0:1:steps),exact_errors/normalization,'linewidth',linesize,'color','black')
hold on
plot((0:1:steps),(min(proj_errors,proj_errors2)+delta_errors)/normalization,'linewidth',linesize,'color',1.2*colors(1,:))  
box on
ax=gca; ax.FontSize=axissize; axis tight
title('Lorenz system','fontsize',fontsize,'interpreter','latex')
xlabel('Time steps','interpreter','latex','fontsize',fontsize)
ylabel('RKHS Error','interpreter','latex','fontsize',fontsize)
exportgraphics(gcf,'lorenz_oscillator_kernel_error_bound.pdf','ContentType','vector','BackgroundColor','none')
