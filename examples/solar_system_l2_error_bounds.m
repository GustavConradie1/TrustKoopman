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
maxx = 10*365;

%%
means = zeros(d,1);
stds = zeros(d,1);
for i=1:d
    means(i)=mean(DATA(i,1:maxx));
    stds(i)=std(DATA(i,1:maxx));
    DATA(i,:)=(DATA(i,:)-means(i))/stds(i);
end

%%
x=DATA(:,1:(maxx-1));
y=DATA(:,2:maxx);
M = maxx-1;

%%
N=200;
cent=cent.';
s=0.0522;
v=4.73;
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

%% svd
r=50;
tic
[U,Sig]=eigs(G,r,'largestabs');
toc
Sig = sqrt(Sig); S2 = diag(1./diag(Sig));
A2 = S2*(U'*A*U)*S2; 
L2 = S2*(U'*L*U)*S2; L2 = (L2+L2')/2;
G2 = eye(r);
K2 = A2;

%% PX_start
x0=DATA(:,maxx);
PX_start=zeros(1,N);
for j=1:N
    r1=s*norm(x0-cent(:,j));
    PX_start(j)=r1^v*besselk(v,r1);
end

%%
NN=1000;

%% evals of kernel
Lambda = exp(-(0:1:(r-1))/(1000)).';

%% get means and variance of norms
K_temp=eye(length(G2));
K_op_avg=zeros(steps+1,1); K_op_avg(1)=1;
pf = parfor_progress(steps);
pfcleanup = onCleanup(@() delete(pf));
for i=1:steps
    K_temp=K2*K_temp;
    KGK=K_temp'*G2*K_temp;
    for j=1:NN
        v=randn(r,1);
        v=v.*Lambda;
        K_op_avg(i+1) = K_op_avg(i+1) + sqrt((v'*KGK*v)/(v'*v));
    end
    K_op_avg(i+1)=K_op_avg(i+1)/NN;
    parfor_progress(pf);
end

%%
exact_errors=zeros(steps+1,36);
error_bounds=zeros(steps+1,36);
exp_errors=zeros(steps+1,36);

%%
pf=parfor_progress(36);
pfcleanup=onCleanup(@() delete(pf));

for kk=1:36
    %%
    coefs_svd = (PX\x(kk,:).');
    coefs = Sig*U'*coefs_svd;
    
    fut_coefs=zeros(r,steps+1);
    fut_coefs(:,1)=coefs;
    for i=1:steps
        fut_coefs(:,i+1)=K2*fut_coefs(:,i);
    end
    
    %% exact errors
    exact_errors(:,kk) = zeros(steps+1,1);
    fut_coefs_svd=U*S2*fut_coefs;
    for i=0:steps
        x_step=DATA(kk,i+1:maxx+i-1);
        exact_errors(i+1,kk)=sqrt(x_step*x_step.'/M-(2/M)*((PX*fut_coefs_svd(:,i+1)).'*x_step.')+fut_coefs(:,i+1)'*fut_coefs(:,i+1));
    end
    
    %% strict error bounds
    [proj_errors,K_op]=kmd_error_bound_first_order_gelfand2(coefs,steps,G2,A2,L2,K2);
    proj_errors2=kmd_error_bound_first_order_gelfand2(coefs,steps,G2,A2,L2,K2);
    disp(proj_errors<=proj_errors2)
    
    %%
    delta = sqrt(x(kk,:)*x(kk,:).'/M-(2/M)*((PX*(U*S2*Sig*U.')*coefs_svd).'*x(kk,:).')+coefs'*coefs);
    delta_errors=delta*K_op;

    error_bounds(:,kk) = delta_errors+proj_errors;

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
    delta_errors_avg=delta*K_op_avg;
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
            exp_errors(i+1,kk)=exp_errors(i+1,kk)+norm(sum(v,2));
        end
        exp_errors(i+1,kk)=exp_errors(i+1,kk)/NN; 
    end
    exp_errors(1,kk)=delta;
    parfor_progress(pf);
end

%%
exact_error=mean(exact_errors,2);
error_bound=mean(error_bounds,2);
exp_error=mean(exp_errors,2);

%% Plot relative forecast errors compared to test data w averaging and std 
figure
semilogy(0:1:steps,exact_error,'linewidth',2)
hold on
plot(0:1:steps,error_bound,'linewidth',2)
plot(0:1:steps,exp_error,'linewidth',2)
box on
ax=gca; ax.FontSize=18; axis tight
title('Pluto-Charon','fontsize',18,'interpreter','latex')
xlabel('Time (days)','interpreter','latex','fontsize',18)
ylabel('$L^2$ norm forecast error','interpreter','latex','fontsize',18)
legend({'Exact','Bound','Expected'},'interpreter','latex','fontsize',16,'location','best')

%%
save('solar_system_results_l2.mat','steps','exact_error','error_bound','exp_error')