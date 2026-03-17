clear

%%%%%%%%%%%% L2 error bounds %%%%%%%%%%%%%
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
N=500;

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
cent(2,:)=50*cent(2,:);
cent(3,:)=50*cent(3,:)+25;
scale=1;
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

%% svd
r=300;
tic
[U,Sig]=eigs(G,r,'largestabs');
toc
Sig = sqrt(Sig); S2 = diag(1./diag(Sig));
A2 = S2*(U'*A*U)*S2; 
L2 = S2*(U'*L*U)*S2; L2 = (L2+L2')/2;
G2 = eye(r);
K2 = A2;

%% coefs for kmd
idx=1;
coefs = (PX\x(:,1:M).'); 
coefs = Sig*U'*coefs(:,idx); 

%% exact evolution of future coefs
steps=20;
fut_coefs=zeros(r,steps+1);
fut_coefs(:,1)=coefs;
for i=1:steps
    fut_coefs(:,i+1)=K2*fut_coefs(:,i); 
end

%% get test data
x0=x_test(:,1);

%% pointwise predictions
PX_test=zeros(1,N);
for j=1:N
    PX_test(j)=rbf(norm(x0-cent(:,j))); 
end
x_pred=PX_test*U*S2*fut_coefs;

%%
fontsize=36;
axissize=28;
linesize=5;
colors=orderedcolors('gem');
figure
plot(dt*(0:1:steps)/0.9056,x_test(idx,1:steps+1),'linewidth',linesize,'color','black')
hold on
plot(dt*(0:1:steps)/0.9056,x_pred,'linewidth',linesize,'color',colors(1,:))  
box on
ax=gca; ax.FontSize=axissize; axis tight
title('Lorenz system','fontsize',fontsize,'interpreter','latex')
xlabel('$Lt$','interpreter','latex','fontsize',fontsize)
ylabel('$x_1$','interpreter','latex','fontsize',fontsize)
legend({'Exact','Predicted'},'interpreter','latex','fontsize',fontsize,'location','best')
exportgraphics(gcf,'big_figure_plots\\lorenz_oscillator_error_bounds.pdf','ContentType','vector','BackgroundColor','none')

%% exact errors
exact_errors = zeros(steps+1,1);
coefs_svd=U*S2*coefs;
fut_coefs_svd=U*S2*fut_coefs;

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
pf=parfor_progress(steps);
pfcleanup=onCleanup(@() delete(pf));
G_a=zeros(N,N,steps);
G_b=zeros(N,N,steps);
for i=1:steps
    for j=1:N
        for l=1:M
            PX_temp(l,j)=rbf(norm(x_test(:,l+i)-cent(:,j)));
        end
    end
    G_a(:,:,i)=PX_temp'*PX2/M;
    G_b(:,:,i)=PX_temp'*PX_temp/M;
    exact_errors(i+1)=sqrt(fut_coefs_svd(:,i+1)'*G_test*fut_coefs_svd(:,i+1)-2*real(fut_coefs_svd(:,i+1)'*G_a(:,:,i).'*coefs_svd)+coefs_svd'*G_b(:,:,i)*coefs_svd);
    parfor_progress(pf);
end

%% strict error bounds
[proj_errors,K_op]=kmd_error_bound_first_order_gelfand2(coefs,steps,G2,A2,L2,K2,1); %enforce unitary
proj_errors2=kmd_error_bound_first_order(coefs,steps,G2,A2,L2,K2,1);

%%
NN=1000;

%% evals of kernel
Lambda = exp(-(0:1:(length(coefs)-1))/(1000)).';

%% get means and variance of norms
K_op_avg=ones(steps+1,1); %enforce unitary

%%
Q=L2-A2'*K2; Q=(Q+Q')/2;
Kg=zeros(length(coefs),steps); Kg(:,1)=coefs;
KgQKg=zeros(steps,1); KgQKg(1)=Kg(:,1)'*Q*Kg(:,1);
for i=2:steps
    Kg(:,i)=K2*Kg(:,i-1);
    KgQKg(i)=Kg(:,i)'*Q*Kg(:,i);
end
KgQKg=real(sqrt(KgQKg));

%%
exp_error2=zeros(steps+1,1);
for i=1:steps
    for j=0:i-1
        exp_error2(i+1)=exp_error2(i+1)+(K_op_avg(2)^j*KgQKg(i-j))^2;
    end
    exp_error2(i+1)=sqrt(exp_error2(i+1));
end

%%
exp_error=zeros(steps+1,1);
pf = parfor_progress(steps);
pfcleanup = onCleanup(@() delete(pf));
for i=1:steps
    temp=zeros(1,i);
    for j=0:i-1
        temp(j+1)=K_op_avg(j+1)*KgQKg(i-j);
    end
    for j=1:NN
        v=randn(length(coefs),i);
        v=v.*Lambda;
        v=v./vecnorm(v);
        v=v.*temp;
        exp_error(i+1)=exp_error(i+1)+norm(sum(v,2));
    end
    exp_error(i+1)=exp_error(i+1)/NN; 
    parfor_progress(pf);
end

%%
normalization=sqrt(coefs'*coefs);

%% Plot relative forecast errors compared to test data w averaging and std 
fontsize=36;
axissize=28;
linesize=5;
colors=orderedcolors('gem');
figure
idx=find(proj_errors>=proj_errors2);
exp_error_combined=exp_error;
exp_error_combined(idx)=exp_error2(idx);
plot((0:1:steps),exact_errors/normalization,'linewidth',linesize,'color','black')
hold on
plot((0:1:steps),min(proj_errors2,proj_errors)/normalization,'linewidth',linesize,'color',colors(1,:))  
plot((0:1:steps),exp_error_combined/normalization,'linewidth',linesize,'color',colors(2,:))
box on
ax=gca; ax.FontSize=axissize; axis tight
title('Lorenz system','fontsize',fontsize,'interpreter','latex')
xlabel('Time steps','interpreter','latex','fontsize',fontsize)
ylabel('$L^2$ Error','interpreter','latex','fontsize',fontsize)
legend({'Exact','Bound','Expected'},'interpreter','latex','fontsize',fontsize,'location','best')

%%
save('lorenz_l2_error_bounds.mat','steps','exact_errors','proj_errors','proj_errors2','exp_error','exp_error2','normalization','exp_error_combined')

%%
clear
addpath(genpath('./algorithms'))
addpath(genpath('./colormaps'))

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
cent(2,:)=50*cent(2,:);
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
steps=20;
PX_temp=zeros(M,N);
pf=parfor_progress(steps);
pfcleanup=onCleanup(@() delete(pf));
G_a=zeros(N,N,steps);
G_b=zeros(N,N,steps);
for i=1:steps
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
K = PX\PY;

%% Plot principal angles and vectors

[theta,U1,U2,J]=KoopAngles(G,A,L);
figure
scatter(1:N,theta,'filled')
title('$\theta$','interpreter','latex','fontsize',18)
ax=gca; ax.FontSize=18;

%% convert U1 to lie entirely inside V
[size1,size2]=size(U1);
V1 = U1(1:size1/2,:) + K*U1(size1/2+1:size1,:);

%% convert U2 to lie entirely inside KV
V2 = K\U2(1:size1/2,:) + U2(size1/2+1:size1,:);

%% 
exact_errors_pad=zeros(1,steps+1);

%% get r smallest angles
r=floor(N/2);
U=V1(:,1:r);
G3=U'*G*U; A3=U'*A*U; L3=U'*L*U; K3=A3;

%% get coefficients of x variable
coefs_pad = (PX*U)\transpose(x(:,1:M));
coefs_pad = coefs_pad(:,1);

%% evolution
fut_coefs_pad=zeros(r,steps+1);
fut_coefs_pad(:,1)=coefs_pad;
for i=1:steps
    fut_coefs_pad(:,i+1)=K3*fut_coefs_pad(:,i);
end

%% exact errors
coefs_pad2=U*coefs_pad;
fut_coefs_pad2=U*fut_coefs_pad;
for i=1:steps
    exact_errors_pad(i+1)=sqrt(fut_coefs_pad2(:,i+1)'*G_test*fut_coefs_pad2(:,i+1)-2*real(fut_coefs_pad2(:,i+1)'*G_a(:,:,i).'*coefs_pad2)+coefs_pad2'*G_b(:,:,i)*coefs_pad2);
end

%% strict error bounds
proj_errors_pad=kmd_error_bound_first_order_gelfand2(coefs_pad,steps,G3,A3,L3,K3);
proj_errors2_pad=kmd_error_bound_first_order(coefs_pad,steps,G3,A3,L3,K3);

%%
normalization2=sqrt(coefs_pad'*coefs_pad);

%% svd
tic
[U,Sig]=eigs(G,r,'largestabs');
toc
Sig = sqrt(Sig); S2 = diag(1./diag(Sig));
A2 = S2*(U'*A*U)*S2; 
L2 = S2*(U'*L*U)*S2; L2 = (L2+L2')/2;
G2 = eye(r);
K2 = A2;

%% coefs for kmd
idx=1;
coefs = (PX\x(:,1:M).'); 
coefs = Sig*U'*coefs(:,idx); 

%% exact evolution of future coefs
fut_coefs=zeros(r,steps+1);
fut_coefs(:,1)=coefs;
for i=1:steps
    fut_coefs(:,i+1)=K2*fut_coefs(:,i);
end

%% 
exact_errors=zeros(1,steps+1);
coefs_svd=U*S2*coefs;
fut_coefs_svd=U*S2*fut_coefs;
for i=1:steps
    exact_errors(i+1)=sqrt(fut_coefs_svd(:,i+1)'*G_test*fut_coefs_svd(:,i+1)-2*real(fut_coefs_svd(:,i+1)'*G_a(:,:,i).'*coefs_svd)+coefs_svd'*G_b(:,:,i)*coefs_svd);
end

%% strict error bounds
proj_errors=kmd_error_bound_first_order_gelfand2(coefs,steps,G2,A2,L2,K2);
proj_errors2=kmd_error_bound_first_order(coefs,steps,G2,A2,L2,K2);

%%
normalization1=sqrt(coefs'*coefs);

%% Plot relative forecast errors compared to test data w averaging and std 
figure
semilogy(0:1:steps,exact_errors/normalization1,'linewidth',2)
hold on
plot(0:1:steps,min(proj_errors,proj_errors2)/normalization1,'linewidth',2)  
plot(0:1:steps,exact_errors_pad/normalization2,'linewidth',2)
plot(0:1:steps,min(proj_errors,proj_errors2_pad)/normalization2,'linewidth',2)  
box on
ax=gca; ax.FontSize=18; axis tight
title('Lorenz oscillator','fontsize',18,'interpreter','latex')
xlabel('Time steps','interpreter','latex','fontsize',18)
ylabel('$L^2$ norm forecast error','interpreter','latex','fontsize',18)
legend({'Exact error','$E_1$','$E_2$','Exact error pad','$E_1$ pad','$E_2$ pad'},'interpreter','latex','fontsize',16,'location','best')

%%
%%%%%%%%%%%% save data %%%%%%%%%%%%%%%%%%
save('big_figure_plots\\lorenz_l2_error_bounds_pad.mat','steps','exact_errors','proj_errors','proj_errors2','exact_errors_pad','proj_errors_pad','proj_errors2_pad','normalization1','normalization2')

%%%%%%%%%% kernel error bounds %%%%%%%%%%%%%
%%
rng(0)
options = odeset('RelTol',1e-10,'AbsTol',1e-10);

%% Set parameters
M=5*10^3;     % number of data points
dt=0.1;    % time step for trajectory sampling
dt2=dt;    % time step for delays
h=dt2/dt;

SIGMA=10;   BETA=8/3;   RHO=28;
ODEFUN=@(t,y) [SIGMA*(y(2)-y(1));y(1).*(RHO-y(3))-y(2);y(1).*y(2)-BETA*y(3)];

TT = 1000; %extra data?

%% Produce the data
tic
Y0=(rand(3,1)-0.5)*4;
[~,x]=ode45(ODEFUN,[0.000001 (1:(100000+(M+h*(TT))))*dt],Y0,options);
x=x(100001:end,:).'; % sample after when on the attractor

Y0=(rand(3,1)-0.5)*4;
[~,x_test]=ode45(ODEFUN,[0.000001 (1:(100000+(M+h*(TT))))*dt],Y0,options);
x_test=x_test(100001:end,:).'; % sample after when on the attractor
toc

%%
scale=0.001/mean(std(x.'));

%% compute operator folding matrices
ker=@(x,t) kernel(x,t,scale);
[G,A,R]=generate_matrices_kernelized(x(:,1:M),x(:,2:M+1),ker);

%% svd
r=M/5;
tic
[U,Sig]=eigs(G,r,'largestabs');
toc
Sig = sqrt(Sig); S2 = diag(1./diag(Sig));
A2 = S2*(U'*A*U)*S2; 
R2 = S2*(U'*R*U)*S2; R2 = (R2+R2')/2;
G2 = eye(r);
K2 = A2;

%% predict future trajectory
x0=x_test(:,1);
Kx0_vals=zeros(M,1);
for i=1:M
    Kx0_vals(i)=ker(x0,x(:,i));
end
coefs=(U*Sig\Kx0_vals);
steps=20;
fut_coefs=zeros(r,steps+1);
fut_coefs(:,1)=coefs;
for i=1:steps
    fut_coefs(:,i+1)=K2*fut_coefs(:,i);
end
x_rkhs=x(:,1:M)*U*S2*conj(fut_coefs);

%%
G_extra=zeros(steps+1,M);
for i=1:steps+1
    for j=1:M
        G_extra(i,j)=ker(x_test(:,i),x(:,j));
    end
end

%%
exact_errors=zeros(steps+1,1);
for i=1:steps+1
    exact_errors(i)=sqrt(ker(x_test(:,i),x_test(:,i))+fut_coefs(:,i)'*fut_coefs(:,i)-2*G_extra(i,:)*U*S2*fut_coefs(:,i));
end

%% error bounds on evolution of K_{x_0}
NN=1000;
delta=sqrt(ker(x0,x0)-coefs'*coefs);

%% strict error bounds
[proj_errors,K_op]=kmd_error_bound_first_order_gelfand2(coefs,steps,G2,A2,R2,K2);
[proj_errors2,K_op2]=kmd_error_bound_first_order(coefs,steps,G2,A2,R2,K2);
delta_errors=delta*K_op;

%%
normalization=sqrt(coefs'*coefs);

%% Plot relative forecast errors compared to test data w averaging and std 
figure
semilogy(0:1:steps,exact_errors/normalization,'linewidth',2)
hold on
plot(0:1:steps,(proj_errors+delta_errors)/normalization,'linewidth',2)  
plot(0:1:steps,(min(proj_errors,proj_errors2)+delta_errors)/normalization,'linewidth',2)  
box on 
ax=gca; ax.FontSize=18; %axis tight
title('Lorenz63 system','fontsize',18,'interpreter','latex')
xlabel('Time steps','interpreter','latex','fontsize',18)
ylabel('RKHS Norm forecast error','interpreter','latex','fontsize',18)
legend({'Exact error','Error bound','Projection error','Delta error'},'interpreter','latex','fontsize',16,'location','best')

%%
Lambda = exp(-(0:1:(length(coefs)-1))/(1000)).';

%% get means and variance of norms
K_temp=eye(length(G2));
K_op_avg=zeros(steps+1,1); K_op_avg(1)=1;
pf = parfor_progress(steps);
pfcleanup = onCleanup(@() delete(pf));
for i=1:steps
    K_temp=K2*K_temp;
    KGK=K_temp'*G2*K_temp;
    for j=1:NN
        v=randn(length(coefs),1);
        v=v.*Lambda;
        K_op_avg(i+1) = K_op_avg(i+1) + sqrt((v'*KGK*v)/(v'*v));
    end
    K_op_avg(i+1)=K_op_avg(i+1)/NN;
    parfor_progress(pf);
end

%%
delta_errors_avg=delta*K_op_avg;

%%
Q=R2-A2'*K2; Q=(Q+Q')/2;
Kg=zeros(length(coefs),steps); Kg(:,1)=coefs;
KgQKg=zeros(steps,1); KgQKg(1)=Kg(:,1)'*Q*Kg(:,1);
for i=2:steps
    Kg(:,i)=K2*Kg(:,i-1);
    KgQKg(i)=Kg(:,i)'*Q*Kg(:,i);
end
KgQKg=real(sqrt(KgQKg));

%%
exp_error=zeros(steps+1,1); exp_error(1)=delta;
pf = parfor_progress(steps);
pfcleanup = onCleanup(@() delete(pf));
for i=1:steps
    temp=zeros(1,i+1);
    temp(i+1)=delta_errors_avg(i+1);
    for j=0:i-1
        temp(j+1)=K_op_avg(j+1)*KgQKg(i-j);
    end
    for j=1:NN
        v=randn(length(coefs),i+1);
        v=v.*Lambda;
        v=v./vecnorm(v);
        v=v.*temp;
        exp_error(i+1)=exp_error(i+1)+norm(sum(v,2));
    end
    exp_error(i+1)=exp_error(i+1)/NN; 
    parfor_progress(pf);
end

%% Plot relative forecast errors compared to test data w averaging and std 
figure
semilogy(0:1:steps,exact_errors/normalization,'linewidth',2)
hold on
plot(0:1:steps,(proj_errors+delta_errors)/normalization,'linewidth',2)  
plot(0:1:steps,(exp_error)/normalization,'linewidth',2)  
box on
ax=gca; ax.FontSize=18;
title('Lorenz63 system','fontsize',18,'interpreter','latex')
xlabel('Time steps','interpreter','latex','fontsize',18)
ylabel('RKHS Norm forecast error','interpreter','latex','fontsize',18)
legend({'Exact error','Error bound','Expected error'},'interpreter','latex','fontsize',16,'location','best')

%%
%%%%%%%%%%%% save data %%%%%%%%%%%%%%%%%%
save('big_figure_plots\\lorenz_kernel_error_bounds.mat','steps','exact_errors','proj_errors','delta_errors','exp_error','delta_errors_avg','normalization')


%% define kernel
function ker=kernel(x,t,scale)
    r=vecnorm(x-t);
    r=scale*r;
    ker=zeros(1,size(x,2));
    ker(r>0)=(r(r>0)).^(1/2).*besselk(-1/2,r(r>0));
    ker(r==0)=sqrt(pi/2);
end