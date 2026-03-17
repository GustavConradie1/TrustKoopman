%%
clear
fontsize=48;
axissize=36;
linesize=5;
colors=orderedcolors('gem');
purple = [150/256 50/256 150/256];

%% generate snapshot data for duffing
rng(0)
time_step=0.05;
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
    plot(x(1,:,i),x(2,:,i),'LineWidth',2,'Color',map(i,:));
end
box on
ax=gca; axis([-2.5 2.5 -2.5 2.5])
ax.FontSize=axissize;
axis off
exportgraphics(gcf,'duffing_attractor.png','ContentType','vector','BackgroundColor','none')

%% principal observables
load('big_figure_plots\\duffing_princ_obs.mat')
figure
vals=real(V1(:,10));
scatter(x(1,:),x(2,:),200,vals,'.','LineWidth',1);
box on
axis([-2.5 2.5 -2.5 2.5])
axis off
clim([mean(vals)-std(vals), mean(vals)+std(vals)]) 
colormap('inferno')
exportgraphics(gcf,'duffing_princobs.png')

%%
load('duffing_l2_error_bounds.mat')
normalization1=sqrt(normalization1);
normalization2=sqrt(normalization2);

%% Plot relative forecast errors compared to test data w averaging and std 
figure
semilogy(0:1:steps,exact_errors/normalization1,'linewidth',linesize,'color','black')
hold on
plot(0:1:steps,min(proj_errors2,proj_errors)/normalization1,'linewidth',linesize,'color',colors(1,:))  
plot(0:1:steps,exp_error2/normalization1,'linewidth',linesize,'color',colors(2,:))  
box on
ax=gca; ax.FontSize=axissize; axis([1 20 0.05 1.1]);
xlabel('Time steps','interpreter','latex','fontsize',fontsize)
ylabel('$L^2$ Error','interpreter','latex','fontsize',fontsize)
legend({'Exact','Bound','Expected$\;\;$'},'interpreter','latex','fontsize',fontsize,'location','best')
exportgraphics(gcf,'duffing_errors.png')

%% Set parameters
M=5*10^3;
dt=0.01;

SIGMA=10;   BETA=8/3;   RHO=28;
ODEFUN=@(t,y) [SIGMA*(y(2)-y(1));y(1).*(RHO-y(3))-y(2);y(1).*y(2)-BETA*y(3)];

%% Produce the data
options = odeset('RelTol',1e-10,'AbsTol',1e-10);
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
exportgraphics(gcf,'lorenz_attractor.png','ContentType','vector','BackgroundColor','none')

%%
load('big_figure_plots\\lorenz_princ_obs.mat')
figure
ax=gca;
vals=real(V1(:,10));
scatter3(ax,x(1,1:M2),x(2,1:M2),x(3,1:M2),200,vals(1:M2),'.');
hold on
box on
grid off
axis([-25 25 -40 40 0 60])
axis off
clim([mean(vals)-2*std(vals), mean(vals)+2*std(vals)])
colormap('inferno')
view([34.5 6.666])
exportgraphics(gcf,'lorenz_princobs.png')

%%
load('lorenz_l2_error_bounds_pad.mat')

%% Plot relative forecast errors compared to test data w averaging and std 
figure
semilogy((0:1:steps),exact_errors/normalization1,'linewidth',linesize,'color','black')
hold on
plot((0:1:steps),min(proj_errors2,proj_errors)/normalization1,'linewidth',linesize,'color',1.2*colors(1,:))
plot((0:1:steps),min(proj_errors2,proj_errors)/normalization1,'linewidth',linesize,'color',colors(2,:))
box on
yticks([10^(-2) 10^(-1)])
xlim tight
ax=gca; ax.FontSize=axissize;
xlabel('Time steps','interpreter','latex','fontsize',fontsize)
ylabel('$L^2$ Error','interpreter','latex','fontsize',fontsize)
legend({'Exact','Bound','Expected$\;$'},'interpreter','latex','fontsize',fontsize,'location','southeast')
exportgraphics(gcf,'lorenz_oscillator_error_avg.png')
    
%% pluto snapshots
load('states_pluto2.mat')
load('states_charon2.mat')
load('states_hydra2.mat')
load('states_styx2.mat')
load('states_nix2.mat')
load('states_kerberos2.mat')

DATA=[states_pluto2 states_charon2 states_hydra2 states_styx2 states_nix2 states_kerberos2];
DATA=DATA.';

%%
idx=1:1000;
figure
scatter3(DATA(7,idx)-DATA(1,idx),DATA(8,idx)-DATA(2,idx),DATA(9,idx)-DATA(3,idx),200,'.')
hold on
scatter3(DATA(13,idx)-DATA(1,idx),DATA(14,idx)-DATA(2,idx),DATA(15,idx)-DATA(3,idx),200,'.')
scatter3(DATA(25,idx)-DATA(1,idx),DATA(26,idx)-DATA(2,idx),DATA(27,idx)-DATA(3,idx),200,'.')
scatter3(DATA(31,idx)-DATA(1,idx),DATA(32,idx)-DATA(2,idx),DATA(33,idx)-DATA(3,idx),200,'.')
scatter3(DATA(19,idx)-DATA(1,idx),DATA(20,idx)-DATA(2,idx),DATA(21,idx)-DATA(3,idx),200,'.')
grid off
box off
axis off
view([-4.1857 21.8417])
[~,objh] = legend({'Charon','Hydra','Nix','Kerberos','Styx'},'interpreter','latex','fontsize',axissize,'location','best');
objhl = findobj(objh, 'type', 'patch'); 
set(objhl, 'Markersize', 20); 
exportgraphics(gcf,'pluto_snapshots.png')

%% principal obs
load('pluto.mat')
figure
idx=1:1000;
ax=gca; 
vals=real(V1(:,3));
scatter3(ax,x(7,idx)*stds(7)+means(7)-DATA(1,idx),x(8,idx)*stds(8)+means(8)-DATA(2,idx),x(9,idx)*stds(9)+means(9)-DATA(3,idx),200,vals(idx),'.','LineWidth',1);
hold on
scatter3(ax,x(13,idx)*stds(13)+means(13)-DATA(1,idx),x(14,idx)*stds(14)+means(14)-DATA(2,idx),x(15,idx)*stds(15)+means(15)-DATA(3,idx),200,vals(idx),'.','LineWidth',1);
scatter3(ax,x(25,idx)*stds(25)+means(25)-DATA(1,idx),x(26,idx)*stds(26)+means(26)-DATA(2,idx),x(27,idx)*stds(27)+means(27)-DATA(3,idx),200,vals(idx),'.','LineWidth',1);
scatter3(ax,x(31,idx)*stds(31)+means(31)-DATA(1,idx),x(32,idx)*stds(32)+means(32)-DATA(2,idx),x(33,idx)*stds(33)+means(33)-DATA(3,idx),200,vals(idx),'.','LineWidth',1);
scatter3(ax,x(19,idx)*stds(19)+means(19)-DATA(1,idx),x(20,idx)*stds(20)+means(20)-DATA(2,idx),x(21,idx)*stds(21)+means(21)-DATA(3,idx),200,vals(idx),'.','LineWidth',1);
clim([mean(vals(idx))-1.5*std(vals(idx)), mean(vals(idx))+1.5*std(vals(idx))]) 
colormap('inferno')
box off
grid off
axis off
view([-4.1857 21.8417])
exportgraphics(gcf,'pluto_princobs.png')

%% pluto error bounds expected errors
load('solar_system_results_l2.mat')
figure
semilogy(0:1:steps,exact_error,'linewidth',linesize,'color','black')
hold on
plot(0:1:steps,error_bound,'linewidth',linesize,'color',colors(1,:))  
plot(0:1:steps,exp_error,'linewidth',linesize,'color',colors(2,:))
box on
ax=gca; ax.FontSize=axissize; axis tight
xlabel('Time (days)','interpreter','latex','fontsize',fontsize)
ylabel('$L^2$ Error','interpreter','latex','fontsize',fontsize)
legend({'Exact','Bound','Expected$\;$'},'interpreter','latex','fontsize',fontsize,'location','best')
exportgraphics(gcf,'pluto_errors.png')

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
exportgraphics(gcf,'cavity_data.png')

%% princobs
idx=5;
u=modes(idx,:);
w = u.*abs(u).^(-0.5);
[~,I2]=sort(real(w),'ascend');
x2=xpt(I2);
I3=I2(x2>-0.7);
y2=ypt(I3);
I4=I3(y2<0.2);
figure
scatter(xpt(I4),ypt(I4),3+2*(xpt(I4)>1),real(w(I4)),'filled');
set(gca,'xticklabel',{[]})
set(gca,'yticklabel',{[]})
s = std( real(w( (xpt(:)<1.5)&(abs(u(:))>0) ) ));
m = mean(real(w( (xpt(:)<1.5)&(abs(u(:))>0)  )));
clim([m-3*s,m+3*s])
colormap('inferno')
ax=gca; ax.FontSize=axissize; axis equal;
grid off
axis off
exportgraphics(gcf,'cavity_princobs.png')

%%
load('cavity_flow_results.mat')
normalization=sqrt(normalization);
figure
semilogy(T(2:steps),exact_error(2:steps)/normalization,'linewidth',linesize,'color','black');
hold on
semilogy(T(2:steps),min(error_bounds(2:steps),error_bounds2(2:steps))/normalization,'linewidth',linesize,'color',1.2*colors(1,:));
semilogy(T(2:steps),exp_error_combined(2:steps)/normalization,'linewidth',linesize,'color',colors(2,:));
legend({'Exact','Bound','Expected$\;$'},'interpreter','latex','fontsize',fontsize,'location','southeast')
xlabel('Time(s)','fontsize',fontsize,'interpreter','latex')
ylabel('$L^2$ Error','fontsize',fontsize,'interpreter','latex')
ax = gca; ax.FontSize = axissize; axis([0.01 1 4*10^(-4) 0.14])
yticks([10^(-3) 10^(-2) 10^(-1)])
box on
exportgraphics(gcf,'cavity_errors.png','ContentType','vector','BackgroundColor','none')

%% define duffing oscillator
function dxdt = duffing(~,x)
    dxdt=zeros(2,1);
    dxdt(1)=x(2);
    dxdt(2)=-0.05*x(2)+x(1)-x(1).^3;
end
