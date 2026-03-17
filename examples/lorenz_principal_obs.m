clear

%%%%%%%%%%%%%%%%%%% get EDMD matrices %%%%%%%%%%%%%%%%  

%% Set parameters
M=5*10^4;     % number of data points
dt=0.1;    % time step for trajectory sampling
dt2=dt;    % time step for delays
h=dt2/dt;

SIGMA=10;   BETA=8/3;   RHO=28;
ODEFUN=@(t,y) [SIGMA*(y(2)-y(1));y(1).*(RHO-y(3))-y(2);y(1).*y(2)-BETA*y(3)];
options = odeset('RelTol',1e-10,'AbsTol',1e-10);

TT = 40;
N=100;

%% Produce the data
tic
Y0=(rand(3,1)-0.5)*4;
[~,x]=ode45(ODEFUN,[0.000001 (1:(100000+(M+h*(TT+N))))*dt],Y0,options);
x=x(100001:end,:).'; % sample after when on the attractor

Y0=(rand(3,1)-0.5)*4;
[~,x_test]=ode45(ODEFUN,[0.000001 (1:(100000+(M+h*(TT))))*dt],Y0,options);
x_test=x_test(100001:end,:).'; % sample after when on the attractor
toc

%%
cent=(rand(3,N)-0.5); 
cent(1,:)=50*cent(1,:);
cent(2,:)=50*cent(2,:); %scale centers to match roughly size of data
cent(3,:)=50*cent(3,:)+25;
scale=0.1;
rbf=@(r) exp(-scale*r);
PX=zeros(M,N);
PY=zeros(M,N);
pf=parfor_progress(N);
pfcleanup=onCleanup(@() delete(pf));
for j=1:N
    for i=1:M
        PX(i,j)=rbf(norm(x(:,i)-cent(:,j)));
        PY(i,j)=rbf(norm(x(:,i+1)-cent(:,j)));
    end
    parfor_progress(pf);
end

%% EDMD

G = (PX'*PX)/M;
A = (PX'*PY)/M;
L = (PY'*PY)/M;

K=PX\PY;

%%%%%%%%%%%%%%%%%%% principal angles and observables %%%%%%%%%%%%%%%

%% Plot principal angles and vectors

[theta,U1,U2,J]=KoopAngles(G,A,L);
[size1,size2]=size(U1);
V1 = U1(1:size1/2,:) + K*U1(size1/2+1:size1,:);
V1 = PX*V1;

%%
M2=10^4;
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
exportgraphics(gcf,'big_figure_plots\\lorenz_princobs.png')

%%%%%%%%%%%%% save values %%%%%%%%%%%%%
%%
save('big_figure_plots\\lorenz_princ_obs.mat','theta','V1','x','M2')
