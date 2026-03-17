clear
%needs chebfun to run

%%%%%%%%%%%%%%%%%%% get EDMD matrices %%%%%%%%%%%%%%%%

%% generate snapshot data
rng(0)
time_step=0.25;
steps=20;
test_steps=20;
no_traj=400;
M=no_traj*steps;
init=4*rand(2,no_traj+1)-2;
x=zeros(2,M);
y=zeros(2,M);
for i=1:no_traj
    [~, x1]=ode45(@(t,x) duffing(t,x),0:time_step:(steps*time_step),init(:,i));
    x1=x1';
    for j=1:steps
        x(:,(i-1)*steps+j)=x1(:,j);
        y(:,(i-1)*steps+j)=x1(:,j+1);
    end
end

%%
N=100;
PX=zeros(M,N);
PY=zeros(M,N);
pf=parfor_progress(N);
pfcleanup=onCleanup(@() delete(pf));
for j=1:sqrt(N)
    for k=1:sqrt(N)
        f1=chebpoly(j-1,[-2.5 2.5]);
        f2=chebpoly(k-1,[-2.5 2.5]);
        for i=1:M
            PX(i,(j-1)*sqrt(N)+k)=f1(x(1,i))*f2(x(2,i));
            PY(i,(j-1)*sqrt(N)+k)=f1(y(1,i))*f2(y(2,i));
        end
        parfor_progress(pf);
    end
end

%% EDMD

G = (PX'*PX)/M;
A = (PX'*PY)/M;
L = (PY'*PY)/M;

K=PX\PY;
%K=PX(1:M,:)\PY(1:M,:);

%%%%%%%%%%%%%%%%%%%%%% principal angles and observables %%%%%%%%%%%%%%%%%%

%% angles
[theta,U1,U2,J]=KoopAngles(G,A,L);
[size1,size2]=size(U1);
V1 = U1(1:size1/2,:) + K*U1(size1/2+1:size1,:);
V1 = PX*V1;

%%
figure
vals=real(V1(:,10));
scatter(x(1,:),x(2,:),200,vals,'.','LineWidth',1);
box on
axis([-2.5 2.5 -2.5 2.5])
axis off
clim([mean(vals)-std(vals), mean(vals)+std(vals)]) 
colormap('inferno')


%%%%%%%%%%%%%%% saves values %%%%%%%%%%%%%%%%%

%%
save('big_figure_plots\\duffing_princ_obs.mat','x','theta','V1')

%% define duffing oscillator
function dxdt = duffing(~,x)
    dxdt=zeros(2,1);
    dxdt(1)=x(2);
    dxdt(2)=-0.05*x(2)+x(1)-x(1).^3;
end
