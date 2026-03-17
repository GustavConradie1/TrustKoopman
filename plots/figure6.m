%%
fontsize=36;
axissize=28;
linesize=5;
colors=orderedcolors('gem');
purple = [170/256 40/256 130/256];

%% Plot errors
figure
p2=plot(0:1:10,[0 0 0 0 1 1 1 1 1 1 1],'linewidth',linesize,'color',1.2*colors(1,:));
hold on
p3=plot(0:1:10,sqrt([0 0 0 0 1 2 3 4 5 6 7]),'linewidth',linesize,'color','magenta');
p1=plot(0:1:10,[0 0 0 0 1 1 1 1 1 1 1],'.','markersize',40,'color','black');
box on
ax=gca; ax.FontSize=axissize;
legend([p1 p2 p3],{'Exact','First-Order Bound$\;\;\;$','Full-Order Bound'},'interpreter','latex','fontsize',fontsize,'location','best')
xlabel('Time steps','interpreter','latex','fontsize',fontsize)
ylabel('$L^2$ Error','interpreter','latex','fontsize',fontsize)
title('System with Lebesgue spectrum','interpreter','latex','fontsize',fontsize)
exportgraphics(gcf,'sharpness_of_inequality.pdf','ContentType','vector','BackgroundColor','none')
