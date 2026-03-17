clear
options = odeset('RelTol',1e-10,'AbsTol',1e-10);

rng(0)
axissize=36;

%% Set parameters
M=5*10^3;     % number of data points
dt=0.01;    % time step for trajectory sampling

SIGMA=10;   BETA=8/3;   RHO=28;
ODEFUN=@(t,y) [SIGMA*(y(2)-y(1));y(1).*(RHO-y(3))-y(2);y(1).*y(2)-BETA*y(3)];

%% Produce the data
tic
Y0=(rand(3,1)-0.5)*4;
[~,x]=ode45(ODEFUN,[0.000001 (1:(100000+M))*dt],Y0,options);
x=x(100001:end,:).'; % sample after when on the attractor

Y0=(rand(3,1)-0.5)*4;
[~,x_test]=ode45(ODEFUN,[0.000001 (1:(100000+(M)))*dt],Y0,options);
x_test=x_test(100001:end,:).'; % sample after when on the attractor
toc

%% plot lorenz attractors
figure
map = colormap(coolwarm(100));
ax=gca;
plot3(ax,x(1,:),x(2,:),x(3,:),'LineWidth',2,'color',map(15,:));
hold on
plot3(ax,x_test(1,:),x_test(2,:),x_test(3,:),'LineWidth',2,'color',map(85,:));
box on
ax.FontSize=axissize; 
axis([-20 20 -30 30 0 50])
axis off
view([34.5 6.666])
exportgraphics(gcf,'lorenz_attractor.pdf','ContentType','vector','BackgroundColor','none')

%%
fontsize=48;
linesize=5;

load('lorenz_pad_vs_svd.mat')
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
exportgraphics(gcf,'lorenz_full_pad_svd.pdf','ContentType','vector','BackgroundColor','none')

%% plot errors
figure
hold on
set(gca, 'YScale', 'log')
plot(r_range,error_svd(:,1),'linewidth',linesize,'color',1.2*colors(1,:))    
hold on
plot(r_range,error_pad(:,1),'linewidth',linesize,'color',colors(2,:))
plot([26 26],[0.22 0.004],'--','linewidth',linesize,'color','black')
box on
ax=gca; ax.FontSize=axissize; axis tight
yticks([10^(-2) 10^(-1)])
xlabel('Number of modes','interpreter','latex','fontsize',fontsize)
ylabel('$L^2$ Error','interpreter','latex','fontsize',fontsize)
title('One-step predictions','interpreter','latex','fontsize',fontsize)
legend({'SVD','PAD'},'interpreter','latex','fontsize',axissize,'location','best')
exportgraphics(gcf,'lorenz_pad_vs_svd_modes.pdf','ContentType','vector','BackgroundColor','none')