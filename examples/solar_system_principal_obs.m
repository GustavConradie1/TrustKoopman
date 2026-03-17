clear
rng(0)

%% load data sets
load('states_pluto2.mat')
load('states_charon2.mat')
load('states_hydra2.mat')
load('states_styx2.mat')
load('states_nix2.mat')
load('states_kerberos2.mat')

DATA=[states_pluto2 states_charon2 states_hydra2 states_styx2 states_nix2 states_kerberos2];
DATA=DATA.';
[d,~]=size(DATA);
no_bodies=d/6;
delta_t=1;
maxx = 40*365;

%%%%%%%%%%%%%%%%%%% plot things %%%%%%%%%%%%%%%%%%%%%

%%
idx=1:1000;
figure
scatter3(DATA(7,idx)-DATA(1,idx),DATA(8,idx)-DATA(2,idx),DATA(9,idx)-DATA(3,idx),100,'.','LineWidth',1)
hold on
scatter3(DATA(13,idx)-DATA(1,idx),DATA(14,idx)-DATA(2,idx),DATA(15,idx)-DATA(3,idx),100,'.','LineWidth',1)
scatter3(DATA(25,idx)-DATA(1,idx),DATA(26,idx)-DATA(2,idx),DATA(27,idx)-DATA(3,idx),100,'.','LineWidth',1)
scatter3(DATA(31,idx)-DATA(1,idx),DATA(32,idx)-DATA(2,idx),DATA(33,idx)-DATA(3,idx),190,'.','LineWidth',1)
scatter3(DATA(19,idx)-DATA(1,idx),DATA(20,idx)-DATA(2,idx),DATA(21,idx)-DATA(3,idx),100,'.','LineWidth',1)
xlabel('$x$','interpreter','latex','fontsize',18)
ylabel('$y$','interpreter','latex','fontsize',18)
zlabel('$z$','interpreter','latex','fontsize',18)
ax=gca; ax.FontSize=18;
view([-4.1857 21.8417])
legend({'Charon','Hydra','Nix','Kerberos','Styx'},'interpreter','latex','fontsize',16,'location','best')
exportgraphics(gcf,'pluto_snapshots.png')

%%%%%%%%%%%%%%%%%%%% setup data %%%%%%%%%%%%%%%%%%%%%%%%%

%%
means = zeros(d,1);
stds = zeros(d,1);
for i=1:d
    means(i)=mean(DATA(i,1:maxx));
    stds(i)=std(DATA(i,1:maxx));
    DATA(i,:)=(DATA(i,:)-means(i))/stds(i);%save these
end

%%
x=DATA(:,1:(maxx-1));
y=DATA(:,2:maxx);
M = maxx-1;

%%
N=200;
cent=cent.';
s=1/5;
v=1;
PX=zeros(M,N);
PY=zeros(M,N);
pf=parfor_progress(N);
pfcleanup=onCleanup(@() delete(pf));
for j=1:N
    for i=1:M
        rx=s*norm(x(:,i)-cent(:,j));
        ry=s*norm(y(:,i)-cent(:,j));
        PX(i,j)=rx^v*besselk(v,rx);
        PY(i,j)=ry^v*besselk(v,ry);
    end
    parfor_progress(pf);
end

%%
steps=20;

%% EDMD

G = (PX'*PX)/M;
A = (PX'*PY)/M;
L = (PY'*PY)/M;
K = pinv(G)*A;

%%%%%%%%%%%%%%%%% principal angles and vectors %%%%%%%%%%%%%%

%% Plot principal angles and vectors

[theta,U1,U2,J]=KoopAngles(G,A,L);
[size1,size2]=size(U1);
V1 = U1(1:size1/2,:) + K*U1(size1/2+1:size1,:);
V1 = PX*V1;

%%
figure
scatter(1:length(theta),theta,'filled')
title('$\theta$','interpreter','latex','fontsize',18)

%% princobs
figure
ax=gca;
vals=real(V1(:,5));
scatter3(ax,x(7,:)*stds(7)+means(7),x(8,:)*stds(8)+means(8),x(9,:)*stds(9)+means(9),60,vals,'.','LineWidth',1);
hold on
scatter3(ax,x(13,:)*stds(13)+means(13),x(14,:)*stds(14)+means(14),x(15,:)*stds(15)+means(15),60,vals,'.','LineWidth',1);
scatter3(ax,x(25,:)*stds(25)+means(25),x(26,:)*stds(26)+means(26),x(27,:)*stds(27)+means(27),60,vals,'.','LineWidth',1);
scatter3(ax,x(31,:)*stds(31)+means(31),x(32,:)*stds(32)+means(32),x(33,:)*stds(33)+means(33),60,vals,'.','LineWidth',1);
scatter3(ax,x(19,:)*stds(19)+means(19),x(20,:)*stds(20)+means(20),x(21,:)*stds(21)+means(21),60,vals,'.','LineWidth',1);
box on
grid on
colormap(turbo); 
legend({'Charon','Hydra','Nix','Kerberos','Styx'},'interpreter','latex','fontsize',16,'location','best')
xlabel('$x$','interpreter','latex','fontsize',18)
ylabel('$y$','interpreter','latex','fontsize',18)
zlabel('$z$','interpreter','latex','fontsize',18)
ax.FontSize=18; 
view([379.5 30.666])
exportgraphics(gcf,'pluto_princobs.png')

%%
save('pluto.mat','x','stds','means','V1')
