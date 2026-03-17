%%
clear
fontsize=36;
axissize=28;
linesize=5;
legendsize=24;

load('solar_system_var_results.mat')
colors=orderedcolors('gem');

%% Plot relative forecast errors compared to test data w averaging and std 
figure
semilogy(0:delta_t:delta_t*steps,exact_errors_rkhs/normalization,'linewidth',linesize,'color','black')
hold on
plot(0:delta_t:delta_t*steps,(proj_errors+delta_errors)/normalization,'linewidth',linesize,'color',1.2*colors(1,:))  
plot(0:delta_t:delta_t*steps,exp_error/normalization,'linewidth',linesize,'color',colors(2,:))  
patch([((0:delta_t:delta_t*steps).');flip((0:delta_t:delta_t*steps).')],real([(exp_error_below);flip(exp_error_above)])/normalization,1*colors(2,:),'FaceAlpha',.2)
box on
ax=gca; ax.FontSize=axissize; axis tight
xticks([0 120 240 360 480 600 720])
title('Pluto-Charon system','fontsize',fontsize,'interpreter','latex')
xlabel('Time (days)','interpreter','latex','fontsize',fontsize)
ylabel('RKHS Error','interpreter','latex','fontsize',fontsize)
legend({'Exact','Bound','Expected$\;$'},'interpreter','latex','fontsize',fontsize,'location','best')
exportgraphics(gcf,'pluto_rkhs_errors.pdf','ContentType','vector','BackgroundColor','none')