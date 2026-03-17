clear

%%
load('Re7500.mat')
TT = 100;

N=200;      % number of basis functions used
M1=2*N;     % number of snapshots to compute the basis
M2=5*N;    % number of snapshots used for ResDMD matrices
M3=TT+M2;    % number of snapshots used for ResDMD matrices (plus future predictions for comparing results)

use_DMD=1;

ind1=1:M1;
ind2=1:M3;
ind3=1:M2;
MM2=length(ind3);

%% Apply ResDMD
if use_DMD~=1
    [PX,PY] = kernel_ResDMD(DATA(:,ind1),DATA(:,ind1+1),DATA(:,ind2),DATA(:,ind2+1),'N',N,'Parallel','off','cut_off',0);
else
    [~,S,V]=svd(transpose(DATA(:,ind1))/sqrt(floor(M1)),'econ');
    PX=transpose(DATA(:,ind2))*V(:,1:N)*diag(1./(diag(S(1:N,1:N))));
    PY=transpose(DATA(:,ind2+1))*V(:,1:N)*diag(1./(diag(S(1:N,1:N))));
    clear S V
end

%%
G = (PX(ind3,:)'*PX(ind3,:))/MM2;
A = (PX(ind3,:)'*PY(ind3,:))/MM2;
L = (PY(ind3,:)'*PY(ind3,:))/MM2;

K=G\A;

%% Koopman and coherent mode decompositions
XI=(PX(ind3,:))\transpose(DATA(:,ind3));

%%
Er_DMD = zeros(TT+1,length(DATA));
v2=XI;
VV=XI;
pf = parfor_progress(TT);
pfcleanup = onCleanup(@() delete(pf));
for jj=1:TT
    v2=K*v2;
    d=PX(ind3,:)*v2-PX((ind3)+jj,:)*VV;
    Er_DMD(jj+1,:)=sqrt(real(dot(d,d))/MM2);
    parfor_progress(pf);
end

%%
tic
idx = 108000;
exact_error=Er_DMD(:,idx);
III=1:144000;

%%
coefs=XI(:,idx);
steps=TT;

%%
[error_bounds,K_op] = kmd_error_bound_first_order_gelfand2(coefs,TT,G,A,L,K);
[error_bounds2,K_op2] = kmd_error_bound_first_order(coefs,TT,G,A,L,K);

%% evals of kernel
Lambda = exp(-(0:1:(length(coefs)-1))/(1000)).';

%% get means and variance of norms
NN=1000;
K_temp=eye(length(G));
K_op_avg=zeros(steps+1,1); K_op_avg(1)=1;
K_op_avg2=zeros(steps+1,1); K_op_avg2(1)=1;
K_op_std=zeros(steps+1,1);
pf = parfor_progress(steps);
pfcleanup = onCleanup(@() delete(pf));
for i=1:steps
    K_temp=K*K_temp;
    KGK=K_temp'*G*K_temp;
    for j=1:NN
        v=randn(length(coefs),1);
        v=v.*(Lambda);
        K_op_avg(i+1) = K_op_avg(i+1) + sqrt((v'*KGK*v)/(v'*G*v));
        K_op_avg2(i+1) = K_op_avg2(i+1) + (v'*KGK*v)/(v'*G*v);
    end
    K_op_avg(i+1)=K_op_avg(i+1)/NN;
    K_op_avg2(i+1)=K_op_avg2(i+1)/NN;
    K_op_std(i+1)=sqrt(K_op_avg2(i+1)-K_op_avg(i+1)^2);
    parfor_progress(pf);
end

%%
Q=L-A'*K; Q=(Q+Q')/2;
Kg=zeros(length(coefs),steps); Kg(:,1)=coefs;
KgQKg=zeros(steps,1); KgQKg(1)=Kg(:,1)'*Q*Kg(:,1);
for i=2:steps
    Kg(:,i)=K*Kg(:,i-1);
    KgQKg(i)=Kg(:,i)'*Q*Kg(:,i);
end
KgQKg=real(sqrt(KgQKg));

%%
k=3;
exp_error2=zeros(steps+1,1);
exp_error2_above=zeros(steps+1,1);
exp_error2_below=zeros(steps+1,1);
for i=1:steps
    for j=0:i-1
        exp_error2(i+1)=exp_error2(i+1)+(K_op_avg(2)^j*KgQKg(i-j))^2;
        exp_error2_above(i+1)=exp_error2_above(i+1)+((K_op_avg(2)+k*K_op_std(2))^j*KgQKg(i-j))^2;
        exp_error2_below(i+1)=exp_error2_below(i+1)+((K_op_avg(2)-k*K_op_std(2))^j*KgQKg(i-j))^2;
    end
    exp_error2(i+1)=sqrt(exp_error2(i+1));
    exp_error2_above(i+1)=sqrt(exp_error2_above(i+1));
    exp_error2_below(i+1)=sqrt(exp_error2_below(i+1));
end

%%
exp_error=zeros(steps+1,1);
exp_error_above=zeros(steps+1,1);
exp_error_above2=zeros(steps+1,1);
exp_error_below=zeros(steps+1,1);
exp_error_below2=zeros(steps+1,1);
pf = parfor_progress(steps);
pfcleanup = onCleanup(@() delete(pf));
for i=1:steps
    temp=zeros(1,i);
    temp_above=zeros(1,i);
    temp_below=zeros(1,i);
    for j=0:i-1
        temp(j+1)=K_op_avg(j+1)*KgQKg(i-j);
        temp_above(j+1)=(K_op_avg(j+1)+k*K_op_std(j+1))*KgQKg(i-j);
        temp_below(j+1)=(K_op_avg(j+1)-k*K_op_std(j+1))*KgQKg(i-j);
    end
    for j=1:NN
        v=randn(length(coefs),i);
        v=v.*Lambda;
        v=v./vecnorm(v);
        v1=v.*temp;
        v2=v.*temp_above;
        v3=v.*temp_below;
        exp_error(i+1)=exp_error(i+1)+norm(sum(v1,2));
        exp_error_above(i+1)=exp_error_above(i+1)+norm(sum(v2,2));
        exp_error_above2(i+1)=exp_error_above2(i+1)+norm(sum(v2,2))^2;
        exp_error_below(i+1)=exp_error_below(i+1)+norm(sum(v3,2));
        exp_error_below2(i+1)=exp_error_below2(i+1)+norm(sum(v3,2))^2;
    end
    parfor_progress(pf);
end
exp_error=exp_error/NN; exp_error_above=exp_error_above/NN; exp_error_above2=exp_error_above2/NN;
exp_error_below=exp_error_below/NN; exp_error_below2=exp_error_below2/NN;
exp_error_above_std=sqrt(exp_error_above2-exp_error_above.^2); exp_error_below_std=sqrt(exp_error_below2-exp_error_below.^2);
exp_error_above=exp_error_above+k*exp_error_above_std; exp_error_below=exp_error_below-k*exp_error_below_std;
toc

%%
steps=length(error_bounds);
normalization=coefs'*G*coefs;

%%
colors=orderedcolors('gem');
exp_error_combined=max(exp_error,exp_error2);
exp_error_above_combined=max(exp_error_above,exp_error2_above);
exp_error_below_combined=min(exp_error_below,exp_error2_below);
figure
semilogy(T(1:steps),exact_error/sqrt(normalization),'linewidth',5,'color','black');
hold on
semilogy(T(1:steps),min(error_bounds,error_bounds2)/sqrt(normalization),'linewidth',5,'color',colors(1,:));
semilogy(T(1:steps),exp_error_combined/sqrt(normalization),'linewidth',5,'color',colors(2,:));
patch([T(2:steps);flip(T(2:steps))],real([exp_error_below_combined(2:end);flip(exp_error_above_combined(2:end))])/sqrt(normalization),colors(2,:),'FaceAlpha',.2)
ax = gca; ax.FontSize = 18;
box on
legend({'Exact','Bound','Expected$\;$'},'interpreter','latex','fontsize',18,'location','northwest')
xlabel('Time(s)','fontsize',24,'interpreter','latex')
ylabel('Error','fontsize',24,'interpreter','latex')
title('Cavity flow','interpreter','latex','fontsize',24)

%%
save('big_figure_plots\\cavity_flow_results.mat','steps','normalization','exact_error','error_bounds','error_bounds2','exp_error_combined','exp_error_above_combined','exp_error_below_combined','T')