clear
fontsize=36;
axissize=28;
linesize=5;
colors=orderedcolors('gem');

%% Plot relative forecast errors compared to test data
load('duffing_l2_error_bounds.mat')
normalization1=sqrt(normalization1);
normalization2=sqrt(normalization2);
figure
semilogy(0:1:steps,exact_errors/normalization1,'linewidth',linesize,'color','black')
hold on
plot(0:1:steps,min(proj_errors2,proj_errors)/normalization1,'linewidth',linesize,'color',1.2*colors(1,:))  
plot(0:1:steps,exact_errors_pad/normalization2,':','linewidth',linesize,'color','black')
plot(0:1:steps,min(proj_errors2_pad,proj_errors_pad)/normalization2,':','linewidth',linesize,'color',1.2*colors(1,:))  
box on
ax=gca; ax.FontSize=axissize; xlim tight
yticks([10^(-2) 10^(-1) 10^0])
title('Duffing oscillator','fontsize',fontsize,'interpreter','latex')
xlabel('Time steps','interpreter','latex','fontsize',fontsize)
ylabel('$L^2$ Error','interpreter','latex','fontsize',fontsize)
legend({'Exact SVD','Bound SVD','Exact PAD $\;$','Bound PAD$\;$'},'interpreter','latex','fontsize',fontsize,'location','southeast')
exportgraphics(gcf,'duffing_oscillator_error_bounds_pad.pdf','ContentType','vector','BackgroundColor','none')

%%
load('lorenz_l2_error_bounds_pad.mat')
figure
semilogy((0:1:steps),exact_errors/normalization1,'linewidth',linesize,'color','black')
hold on
plot((0:1:steps),min(proj_errors2,proj_errors)/normalization1,'linewidth',linesize,'color',1.2*colors(1,:))  
plot((0:1:steps),exact_errors_pad/normalization2,':','linewidth',linesize,'color','black')
plot((0:1:steps),min(proj_errors2_pad,proj_errors_pad)/normalization2,':','linewidth',linesize,'color',1.2*colors(1,:))  
box on
ax=gca; ax.FontSize=axissize; xlim tight
yticks([10^(-2) 10^(-1) 10^0])
title('Lorenz system','fontsize',fontsize,'interpreter','latex')
xlabel('Time steps','interpreter','latex','fontsize',fontsize)
ylabel('$L^2$ Error','interpreter','latex','fontsize',fontsize)
exportgraphics(gcf,'lorenz_oscillator_error_bounds_pad.pdf','ContentType','vector','BackgroundColor','none')

