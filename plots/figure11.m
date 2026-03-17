clear
load('cavity_flow_results.mat')
fontsize=36;
axissize=28;
linesize=5;
colors=orderedcolors('gem');
purple = [150/256 50/256 150/256];
normalization=sqrt(normalization);

%%
figure
semilogy(T(2:steps),exact_error(2:steps)/normalization,'linewidth',linesize,'color','black');
hold on
semilogy(T(2:steps),min(error_bounds(2:steps),error_bounds2(2:steps))/normalization,'linewidth',linesize,'color',1.2*colors(1,:));
semilogy(T(2:steps),exp_error_combined(2:steps)/normalization,'linewidth',linesize,'color',colors(2,:));
patch([T(2:steps);flip(T(2:steps))],real([exp_error_below_combined(2:end);flip(exp_error_above_combined(2:end))])/normalization,colors(2,:),'FaceAlpha',.2)
legend({'Exact','Bound','Expected$\;$'},'interpreter','latex','fontsize',fontsize,'location','southeast')
xlabel('Time(s)','fontsize',fontsize,'interpreter','latex')
ylabel('$L^2$ Error','fontsize',fontsize,'interpreter','latex')
ax = gca; ax.FontSize = axissize;
axis([0.01 1 4*10^(-4) 0.14])
yticks([10^(-3) 10^(-2) 10^(-1)])
box on
title('Cavity flow','interpreter','latex','fontsize',fontsize)
exportgraphics(gcf,'cavity_flow_expected_error.pdf','ContentType','vector','BackgroundColor','none')

