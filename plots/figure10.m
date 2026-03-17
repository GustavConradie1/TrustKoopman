clear
fontsize=36;
axissize=28;
linesize=5;
colors=orderedcolors('gem');

%%
load('duffing_l2_error_bounds.mat')
normalization1=sqrt(normalization1);
normalization2=sqrt(normalization2);

%% Plot relative forecast errors compared to test data
figure
semilogy(0:1:steps,exact_errors/normalization1,'linewidth',linesize,'color','black')
hold on
plot(0:1:steps,min(proj_errors2,proj_errors)/normalization1,'linewidth',linesize,'color',1.2*colors(1,:))  
plot(0:1:steps,exp_error2/normalization1,'linewidth',linesize,'color',colors(2,:))  
box on
ax=gca; ax.FontSize=axissize; axis([1 20 0.05 1.1]);
title('Duffing oscillator','fontsize',fontsize,'interpreter','latex')
xlabel('Time steps','interpreter','latex','fontsize',fontsize)
ylabel('$L^2$ Error','interpreter','latex','fontsize',fontsize)
legend({'Exact','Bound','Expected$\;\;$'},'interpreter','latex','fontsize',fontsize,'location','northwest')
exportgraphics(gcf,'duffing_oscillator_error_avg.pdf','ContentType','vector','BackgroundColor','none')

%%
load('duffing_kernel_error_bounds.mat')
normalization=sqrt(normalization);

%% 
figure
semilogy(0:1:steps,exact_errors/normalization,'linewidth',linesize,'color','black')
hold on
plot(0:1:steps,min(proj_errors2,proj_errors)/normalization,'linewidth',linesize,'color',1.2*colors(1,:))  
plot(0:1:steps,exp_error/normalization,'linewidth',linesize,'color',colors(2,:))  
box on
ax=gca; ax.FontSize=axissize;
title('Duffing oscillator','fontsize',fontsize,'interpreter','latex')
xlabel('Time steps','interpreter','latex','fontsize',fontsize)
ylabel('RKHS Error','interpreter','latex','fontsize',fontsize)
legend({'Exact','Bound','Expected$\;$'},'interpreter','latex','fontsize',fontsize,'location','northwest')
exportgraphics(gcf,'duffing_oscillator_kernel_exp_error.pdf','ContentType','vector','BackgroundColor','none')

%%
load('lorenz_l2_error_bounds_pad.mat')

%% 
figure
semilogy((0:1:steps),exact_errors/normalization1,'linewidth',linesize,'color','black')
hold on
plot((0:1:steps),min(proj_errors2,proj_errors)/normalization1,'linewidth',linesize,'color',1.2*colors(1,:))
plot((0:1:steps),min(proj_errors2,proj_errors)/normalization1,'linewidth',linesize,'color',colors(2,:))
box on
yticks([10^(-2) 10^(-1) 10^0])
xlim tight
ax=gca; ax.FontSize=axissize;
title('Lorenz system','fontsize',fontsize,'interpreter','latex')
xlabel('Time steps','interpreter','latex','fontsize',fontsize)
ylabel('$L^2$ Error','interpreter','latex','fontsize',fontsize)
legend({'Exact','Bound','Expected$\;$'},'interpreter','latex','fontsize',fontsize,'location','southeast')
exportgraphics(gcf,'lorenz_oscillator_error_avg.pdf','ContentType','vector','BackgroundColor','none')

%%
load('lorenz_kernel_error_bounds.mat')

%%
figure
semilogy((0:1:steps),exact_errors/normalization,'linewidth',linesize,'color','black')
hold on
plot((0:1:steps),min(proj_errors,proj_errors2)/normalization,'linewidth',linesize,'color',1.2*colors(1,:))  
plot((0:1:steps),exp_error/normalization,'linewidth',linesize,'color',colors(2,:))  
box on
ax=gca; ax.FontSize=axissize; axis tight
title('Lorenz system','fontsize',fontsize,'interpreter','latex')
xlabel('Time steps','interpreter','latex','fontsize',fontsize)
ylabel('RKHS Error','interpreter','latex','fontsize',fontsize)
exportgraphics(gcf,'lorenz_oscillator_kernel_exp_error.pdf','ContentType','vector','BackgroundColor','none')
