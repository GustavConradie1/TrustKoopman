%%
clear
fontsize=36;
axissize=28;
linesize=5;

load('duffing_pad_vs_svd.mat')
colors=orderedcolors('gem');

%% plot errors
figure
semilogy(0:1:test_steps,exact_errors,'linewidth',linesize,'color','black')
hold on
plot(0:1:test_steps,exact_errors_svd,'linewidth',linesize,'color',1.2*colors(1,:))    
plot(0:1:test_steps,exact_errors_pad,'linewidth',linesize,'color',colors(2,:))    
box on
ax=gca; ax.FontSize=axissize; xlim tight
xlabel('Time steps','interpreter','latex','fontsize',fontsize)
ylabel('$L^2$ Error','interpreter','latex','fontsize',fontsize)
title('Forecast predictions','interpreter','latex','fontsize',fontsize)
legend({'Full','SVD','PAD'},'interpreter','latex','fontsize',axissize,'location','best')
exportgraphics(gcf,'duffing_full_pad_svd.pdf','ContentType','vector','BackgroundColor','none')

%% plot errors
figure
hold on
set(gca, 'YScale', 'log')
plot(r_range,error_svd(:,1),'linewidth',linesize,'color',1.2*colors(1,:))    
plot(r_range,error_pad(:,1),'linewidth',linesize,'color',colors(2,:))
plot([49 49],[10^0 10^(-4)],'--','linewidth',linesize,'color','black')
box on
ax=gca; ax.FontSize=axissize; xlim tight
xlabel('Number of modes','interpreter','latex','fontsize',fontsize)
ylabel('$L^2$ Error','interpreter','latex','fontsize',fontsize)
title('One-step predictions','interpreter','latex','fontsize',fontsize)
legend({'SVD','PAD'},'interpreter','latex','fontsize',axissize,'location','best')
exportgraphics(gcf,'duffing_pad_vs_svd_modes.pdf','ContentType','vector','BackgroundColor','none')
