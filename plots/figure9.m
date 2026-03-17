%% plot cavity snapshots
load('cavity.mat')
TT=5000;
u=snapshot;
w = u.*abs(u).^(-0.5);
[~,I2]=sort(real(w),'ascend');
x2=xpt(I2);
I3=I2(x2>-0.5);
y2=ypt(I3);
I4=I3(y2<0.2);
figure
scatter(xpt(I4),ypt(I4),3+2*(xpt(I4)>1),real(w(I4)),'filled');
set(gca,'xticklabel',{[]})
set(gca,'yticklabel',{[]})
colormap(coolwarm); 
s = std( real(w( (xpt(:)<1.5)&(abs(u(:))>0) ) ));
m = mean(real(w( (xpt(:)<1.5)&(abs(u(:))>0)  )));
clim([m-3*s,m+3*s])
ax=gca; ax.FontSize=axissize; axis equal;
grid off
axis off
exportgraphics(gcf,'cavity_data.pdf')

%%
clear
fontsize=36;
axissize=28;
linesize=5;
colors=orderedcolors('gem');
purple = [150/256 50/256 150/256];
normalization=sqrt(normalization);

%% just true error and error bounds 
figure
purple = [150/256 50/256 150/256];
semilogy(T(2:steps),exact_error(2:steps)/normalization,'linewidth',linesize,'color','black');
hold on
semilogy(T(2:steps),min(error_bounds(2:steps),error_bounds2(2:steps))/normalization,'linewidth',linesize,'color',1.2*colors(1,:));
annotation('doublearrow',[0.35 0.35],[0.44 0.61],'color',purple,'linewidth',7,'head1width',20,'head2width',20,'head1length',20,'head2length',20)
annotation('doublearrow',[0.75 0.75],[0.58 0.8],'color',purple,'linewidth',7,'head1width',20,'head2width',20,'head1length',20,'head2length',20)
text(0.33,1.3*10^(-2),'Pessimistic bound','fontsize',32,'interpreter','latex','color',purple)
ax=gca; ax.FontSize=axissize;
legend({'Exact','Bound$\;\;$'},'interpreter','latex','fontsize',fontsize,'location','southeast')
xlabel('Time(s)','fontsize',fontsize,'interpreter','latex')
ylabel('$L^2$ Error','fontsize',fontsize,'interpreter','latex')
box on
title('Cavity flow','interpreter','latex','fontsize',fontsize)
axis([0.01 1 4*10^(-4) 0.14])
yticks([10^(-3) 10^(-2) 10^(-1)])
exportgraphics(gcf,'cavity_flow_just_bound.pdf','ContentType','vector','BackgroundColor','none')

