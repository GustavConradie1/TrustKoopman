%%
clear
fontsize=48;
axissize=36;
linesize=5;
legendsize=36;

load('solar_system_results.mat')
colors=orderedcolors('gem');

%% plot exact versus predicted trajectories
for i=19:24
    figure
    box on
    hold on
    p1=plot(0:delta_t:delta_t*steps,stds(i)*x_rkhs(i,:)+means(i)-(stds(i-18)*x_rkhs(i-18,:)+means(i-18)),'linewidth',linesize,'color','magenta');
    p2=plot(0:delta_t:delta_t*steps,stds(i)*real_data(i,:)+means(i)-(stds(i-18)*real_data(i-18,:)+means(i-18)),'--','linewidth',linesize,'color','black');
    xlabel('Time (days)','interpreter','latex','fontsize',fontsize)
        xticks([0 120 240 360 480 600 720])
    ax=gca; ax.FontSize=axissize; xlim('tight')
    if mod(i,6)==1
        title('Styx, $x$-coordinate','interpreter','latex','fontsize',fontsize)
        legend([p2 p1],'Exact','Forecast','interpreter','latex','fontsize',legendsize,'location','southeast')
        exportgraphics(gcf,'pluto_real_vs_pred_x.pdf','ContentType','vector','BackgroundColor','none')
    elseif mod(i,6)==2
        title('Styx, $y$-coordinate','interpreter','latex','fontsize',fontsize)
        exportgraphics(gcf,'pluto_real_vs_pred_y.pdf','ContentType','vector','BackgroundColor','none')
    elseif mod(i,6)==3
        title('Styx, $z$-coordinate','interpreter','latex','fontsize',fontsize)
        exportgraphics(gcf,'pluto_real_vs_pred_z.pdf','ContentType','vector','BackgroundColor','none')
    elseif mod(i,6)==4
        title('Styx, $\dot{x}$-coordinate','interpreter','latex','fontsize',fontsize)
        legend([p2 p1],'Exact','Forecast','interpreter','latex','fontsize',legendsize,'location','southeast')
        exportgraphics(gcf,'pluto_real_vs_pred_xdot.pdf','ContentType','vector','BackgroundColor','none')
    elseif mod(i,6)==5
        title('Styx, $\dot{y}$-coordinate','interpreter','latex','fontsize',fontsize)
        exportgraphics(gcf,'pluto_real_vs_pred_ydot.pdf','ContentType','vector','BackgroundColor','none')
    elseif mod(i,6)==0
        title('Styx, $\dot{z}$-coordinate','interpreter','latex','fontsize',fontsize)
        exportgraphics(gcf,'pluto_real_vs_pred_zdot.pdf','ContentType','vector','BackgroundColor','none')
    end
end

%% Plot pointwise errors bounds compared to test data
for i=19:24
    figure
    semilogy(0:delta_t:delta_t*steps,abs(real_data(i,:)-x_rkhs(i,:)),':','linewidth',linesize,'color','black')
    hold on
    plot(0:delta_t:delta_t*steps,g_norms_sup(i)*(proj_errors+delta_errors),':','linewidth',linesize,'color',1.2*colors(1,:))  
    plot(0:delta_t:delta_t*steps,g_norms_exp(i)*(error),':','linewidth',linesize,'color',colors(2,:))  
    box on
    ax=gca; ax.FontSize=axissize; xlim('tight')
    xlabel('Time (days)','interpreter','latex','fontsize',fontsize)
    xticks([0 120 240 360 480 600 720])
    yticks([10^(-5) 10^(-4) 10^(-3) 10^(-2) 10^(-1) 10^0])
    if mod(i,6)==1
        title('Styx, $x$-coordinate','interpreter','latex','fontsize',fontsize)
        legend({'Exact','Bound','Expected$\;\;$'},'interpreter','latex','fontsize',legendsize,'location','southeast');
        exportgraphics(gcf,'pluto_errors_x.pdf','ContentType','vector','BackgroundColor','none')
    elseif mod(i,6)==2
        title('Styx, $y$-coordinate','interpreter','latex','fontsize',fontsize)
        exportgraphics(gcf,'pluto_errors_y.pdf','ContentType','vector','BackgroundColor','none')
    elseif mod(i,6)==3
        title('Styx, $z$-coordinate','interpreter','latex','fontsize',fontsize)
        exportgraphics(gcf,'pluto_errors_z.pdf','ContentType','vector','BackgroundColor','none')
    elseif mod(i,6)==4
        title('Styx, $\dot{x}$-coordinate','interpreter','latex','fontsize',fontsize)
                legend({'Exact','Bound','Expected$\;\;$'},'interpreter','latex','fontsize',legendsize,'location','southeast');
        exportgraphics(gcf,'pluto_errors_xdot.pdf','ContentType','vector','BackgroundColor','none')
    elseif mod(i,6)==5
        title('Styx, $\dot{y}$-coordinate','interpreter','latex','fontsize',fontsize)
        exportgraphics(gcf,'pluto_errors_ydot.pdf','ContentType','vector','BackgroundColor','none')
    elseif mod(i,6)==0
        title('Styx, $\dot{z}$-coordinate','interpreter','latex','fontsize',fontsize)
        exportgraphics(gcf,'pluto_errors_zdot.pdf','ContentType','vector','BackgroundColor','none')
    end
end
