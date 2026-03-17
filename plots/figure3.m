%%
clear
load('duffing_principal_angle.mat')

%%
fontsize=48;
axissize=36;
dotsize=70;
linesize=5;

%% get data to plot attractor
rng(0)
time_step=0.05; %0.25
steps=500;
no_traj=20;
init=4*rand(2,no_traj)-2;
x=zeros(2,steps+1,no_traj);
for i=1:no_traj
    [~, x1]=ode45(@(t,x) duffing(t,x),0:time_step:(steps*time_step),init(:,i));
    x(:,:,i)=x1';
end

%% duffing attractors   
figure
map = colormap(coolwarm(no_traj));
hold on
for i=1:no_traj
    plot(x(1,:,i),x(2,:,i),'LineWidth',1,'Color',map(i,:));
end
box on
ax=gca; axis([-2.5 2.5 -2.5 2.5])
ax.FontSize=axissize;
xlabel('$x$','interpreter','latex','fontsize',fontsize)
ylabel('$y$','interpreter','latex','fontsize',fontsize,'Rotation', 0)
title('Duffing attractor','interpreter','latex','fontsize',fontsize)
exportgraphics(gcf,'duffing_attractor.pdf','ContentType','vector','BackgroundColor','none')

%% Plot principal angles and vectors
N=length(theta);
figure
hold on
scatter(1:N,theta,dotsize,'black','filled')
yticks([0 pi/4 pi/2])
yticklabels({'$\;\;0$','$\;\;\pi/4$','$\;\;\pi/2$'})
set(gca,'TickLabelInterpreter','latex')
box on
ax=gca; ax.FontSize=axissize; axis([0 100 0 pi/2])
xlabel('Mode','interpreter','latex','fontsize',fontsize)
ylabel('$\theta$','interpreter','latex','fontsize',fontsize,'Rotation', 0)
title('Principal angles','interpreter','latex','fontsize',fontsize)
exportgraphics(gcf,'duffing_full_angles.pdf','ContentType','vector','BackgroundColor','none')

%% best piecewise linear 2 components fits of lines
[coef1,bkpt1,~]=fitBogartz(log(theta(2:end)),log(initial_errors(2:end)),0); 
r1=linspace(theta(2),exp(bkpt1(1)),100); 
s1=linspace(exp(bkpt1(1)),theta(end),100);

%% plot
figure
hold on
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
colors=orderedcolors('gem');
scatter(theta(2:end),initial_errors(2:end),dotsize,1.2*colors(1,:),'filled')
plot(r1,exp(coef1(1)*log(r1)+coef1(2)),'linewidth',linesize,'color','black')
plot(s1,exp(coef1(3)*log(s1)+coef1(4)),'linewidth',linesize,'color','black')
box on
ax=gca; ax.FontSize=axissize;
xlabel('$\theta$','interpreter','latex','fontsize',fontsize)
ylabel('$L^2$ Error','interpreter','latex','fontsize',fontsize)
title('Principal observable prediction','interpreter','latex','fontsize',fontsize)
xticks([10^(-4) 10^(-2) 10^(0)])
yticks([10^(-2) 10^(-1) 10^0])
xlim tight
ylim([0.003 1])

exportgraphics(gcf,'duffing_angles_vs_errors.pdf','ContentType','vector','BackgroundColor','none')


%% define duffing oscillator
function dxdt = duffing(~,x)
    dxdt=zeros(2,1);
    dxdt(1)=x(2);
    dxdt(2)=-0.05*x(2)+x(1)-x(1).^3;
end
