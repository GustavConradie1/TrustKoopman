clear

options = odeset('RelTol',1e-10,'AbsTol',1e-10);
rng(0)

%% Set parameters
M=5*10^4;     % number of data points
dt=0.1;    % time step for trajectory sampling
dt2=dt;    % time step for delays
h=dt2/dt;

SIGMA=10;   BETA=8/3;   RHO=28;
ODEFUN=@(t,y) [SIGMA*(y(2)-y(1));y(1).*(RHO-y(3))-y(2);y(1).*y(2)-BETA*y(3)];

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
K=pinv(G)*A;

%% get coefficients of full state observable for full model
coefs = PX\transpose(x(:,1:M));
coefs=coefs(:,1); 
normalization=sqrt(coefs'*G*coefs);

%% evolution
test_steps=20;
fut_coefs=zeros(N,test_steps+1);
fut_coefs(:,1)=coefs;
for i=1:test_steps
    fut_coefs(:,i+1)=K*fut_coefs(:,i);
end

%% exact errors of observables
x_future=zeros(3,M,test_steps+1);
for i=1:M
    [~, x_temp]=ode45(@(t,x) lorenz(t,x),0:dt:test_steps*dt,x(:,i));
    x_future(:,i,:)=x_temp.';
end

%% get another G matrix for testing comparisons w new quadrature
PX2=zeros(M,N);
pf=parfor_progress(N);
pfcleanup=onCleanup(@() delete(pf));
for j=1:N
    for i=1:M
        PX2(i,j)=rbf(norm(x_test(:,i)-cent(:,j)));
    end
    parfor_progress(pf);
end
G_test = (PX2'*PX2)/M;

%% 
PX_temp=zeros(M,N);
pf=parfor_progress(test_steps);
pfcleanup=onCleanup(@() delete(pf));
G_a=zeros(N,N,test_steps);
G_b=zeros(N,N,test_steps);
for i=1:test_steps
    for j=1:N
        for l=1:M
            PX_temp(l,j)=rbf(norm(x_test(:,l+i)-cent(:,j)));
        end
    end
    G_a(:,:,i)=PX_temp'*PX2/M;
    G_b(:,:,i)=PX_temp'*PX_temp/M;
    parfor_progress(pf);
end

%%
exact_errors=zeros(test_steps+1,1);
for i=1:test_steps
    exact_errors(i+1)=sqrt(fut_coefs(:,i+1)'*G_test*fut_coefs(:,i+1)-2*real(fut_coefs(:,i+1)'*G_a(:,:,i).'*coefs)+coefs'*G_b(:,:,i)*coefs);
end
exact_errors=exact_errors/normalization;

%% Plot principal angles and vectors

[theta,U1,U2,J]=KoopAngles(G,A,L);

%% convert U1 to lie entirely inside V
[size1,size2]=size(U1);
V1 = U1(1:size1/2,:) + K*U1(size1/2+1:size1,:);

%% convert U2 to lie entirely inside KV
V2 = K\U2(1:size1/2,:) + U2(size1/2+1:size1,:);

%% truncate via principle angles
idx=find(theta<0.1); r=length(idx);
U=V1(:,idx);
G2=U'*G*U; A2=U'*A*U; L2=U'*L*U; K2=A2;

%% coefs 
coefs = (PX*U)\transpose(x(:,1:M));
coefs = coefs(:,1);
normalization_pad=sqrt(coefs'*coefs);

%% evolution
fut_coefs=zeros(r,test_steps+1);
fut_coefs(:,1)=coefs;
for i=1:test_steps
    fut_coefs(:,i+1)=K2*fut_coefs(:,i);
end

%%
coefs_pad=U*coefs;
fut_coefs_pad=U*fut_coefs;
exact_errors_pad=zeros(test_steps+1,1);
for i=1:test_steps
    exact_errors_pad(i+1)=sqrt(fut_coefs_pad(:,i+1)'*G_test*fut_coefs_pad(:,i+1)-2*real(fut_coefs_pad(:,i+1)'*G_a(:,:,i).'*coefs_pad)+coefs_pad'*G_b(:,:,i)*coefs_pad);
end
exact_errors_pad=exact_errors_pad/normalization_pad;

%% truncate via usual SVD
[U,Sig]=eigs(G,r,'largestabs');
Sig = sqrt(Sig); S2 = diag(1./diag(Sig));
A3 = S2*(U'*A*U)*S2; 
L3 = S2*(U'*L*U)*S2; L3 = (L3+L3')/2;
G3 = eye(r);
K3 = A3;

%% coefs 
coefs = (PX*U*S2)\transpose(x(:,1:M));
coefs = coefs(:,1);
normalization_svd=sqrt(coefs'*coefs);

%% evolution
fut_coefs=zeros(r,test_steps+1);
fut_coefs(:,1)=coefs;
for i=1:test_steps
    fut_coefs(:,i+1)=K3*fut_coefs(:,i);
end

%%
exact_errors_svd=zeros(test_steps+1,1);
coefs_svd=U*S2*coefs;
fut_coefs_svd=U*S2*fut_coefs;
for i=1:test_steps
    exact_errors_svd(i+1)=sqrt(fut_coefs_svd(:,i+1)'*G_test*fut_coefs_svd(:,i+1)-2*real(fut_coefs_svd(:,i+1)'*G_a(:,:,i).'*coefs_pad)+coefs_pad'*G_b(:,:,i)*coefs_pad);
end
exact_errors_svd=exact_errors_svd/normalization_svd;

%% plot errors
figure
semilogy(0:1:test_steps,exact_errors,'linewidth',2)
hold on
plot(0:1:test_steps,exact_errors_pad,'linewidth',2)    
plot(0:1:test_steps,exact_errors_svd,'linewidth',2)    
box on
ax=gca; ax.FontSize=18; axis tight
xlabel('Time steps','interpreter','latex','fontsize',18)
ylabel('Relative Forecast Error','interpreter','latex','fontsize',18)
legend({'Exact error','Exact error PAD','Exact error SVD'},'interpreter','latex','fontsize',16,'location','best')

%% 
r_range=5:1:size(V1,2);
l=length(r_range);
error_pad=zeros(l,test_steps);
error_svd=zeros(l,test_steps);

for ii=1:l
    ii
    %% get r smallest angles
    r=r_range(ii);
    U=V1(:,1:r);
    G2=U'*G*U; A2=U'*A*U; L2=U'*L*U; K2=A2;
    
    %% get coefficients of x variable
    coefs = (PX*U)\transpose(x(:,1:M));
    coefs = coefs(:,1);
    normalization_pad=sqrt(coefs'*coefs);

    %% evolution
    fut_coefs=zeros(r,test_steps+1);
    fut_coefs(:,1)=coefs;
    for i=1:test_steps
        fut_coefs(:,i+1)=K2*fut_coefs(:,i);
    end

    %% exact errors
    coefs_pad=U*coefs;
    fut_coefs_pad=U*fut_coefs;
    for i=1:test_steps
        error_pad(ii,i)=sqrt(fut_coefs_pad(:,i+1)'*G_test*fut_coefs_pad(:,i+1)-2*real(fut_coefs_pad(:,i+1)'*G_a(:,:,i).'*coefs_pad)+coefs_pad'*G_b(:,:,i)*coefs_pad)/normalization_pad;
    end

    %% truncate via usual SVD
    [U,Sig]=eigs(G,r,'largestabs');
    Sig = sqrt(Sig); S2 = diag(1./diag(Sig));
    A3 = S2*(U'*A*U)*S2; 
    L3 = S2*(U'*L*U)*S2; L3 = (L3+L3')/2;
    G3 = eye(r);
    K3 = A3;
    
    %% coefs 
    coefs = (PX*U*S2)\transpose(x(:,1:M));
    coefs = coefs(:,1);
    normalization_svd=sqrt(coefs'*coefs);
    
    %% evolution
    fut_coefs=zeros(r,test_steps+1);
    fut_coefs(:,1)=coefs;
    for i=1:test_steps
        fut_coefs(:,i+1)=K3*fut_coefs(:,i);
    end
        
    %%
    coefs_svd=U*S2*coefs;
    fut_coefs_svd=U*S2*fut_coefs;

    for i=1:test_steps
        error_svd(ii,i)=sqrt(fut_coefs_svd(:,i+1)'*G_test*fut_coefs_svd(:,i+1)-2*real(fut_coefs_svd(:,i+1)'*G_a(:,:,i).'*coefs_pad)+coefs_pad'*G_b(:,:,i)*coefs_pad)/normalization_svd;
    end
end

%% plot errors
figure
hold on
set(gca, 'YScale', 'log')
plot(r_range,error_pad(:,1),'linewidth',2)
hold on
plot(r_range,error_svd(:,1),'linewidth',2)    
box on
ax=gca; ax.FontSize=18; axis tight
xlabel('Number of modes','interpreter','latex','fontsize',18)
ylabel('One step prediction error','interpreter','latex','fontsize',18)
legend({'Error PAD','Error SVD'},'interpreter','latex','fontsize',16,'location','best')

%% save
save('big_figure_plots\\lorenz_pad_vs_svd.mat','test_steps','exact_errors','exact_errors_svd','exact_errors_pad','r_range','error_pad','error_svd')

%% define lorenz system
function dxdt = lorenz(~,x)
    sigma=10; rho=28; beta=8/3;
    dxdt=zeros(3,1);
    dxdt(1)=sigma*(x(2)-x(1));
    dxdt(2)=x(1)*(rho-x(3))-x(2);
    dxdt(3)=x(1)*x(2)-beta*x(3);
end
