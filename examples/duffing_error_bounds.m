clear

%% requires chebfun to be installed

%%%%%%%%%%%%%% L2 error bounds %%%%%%%%%%%%%

rng(0)

%% generate snapshot data
time_step=0.25;
steps=20;
no_traj=100;
M=no_traj*steps;
init=4*(rand(2,no_traj+1)-0.5);
x=zeros(2,M);
y=zeros(2,M);
for i=1:no_traj
    [~, x1]=ode45(@(t,z) duffing(t,z),0:time_step:(steps*time_step),init(:,i));
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
K = pinv(G)*A;

%% svd
r=N/2;
tic
[U,Sig]=eigs(G,r,'largestabs');
toc
Sig = sqrt(Sig); S2 = diag(1./diag(Sig));
A2 = S2*(U'*A*U)*S2; 
L2 = S2*(U'*L*U)*S2; L2 = (L2+L2')/2;
G2 = eye(r);
K2 = A2;

%% coefs for kmd
coefs = (PX\x.');
coefs = Sig*U'*coefs(:,1);

%% get test data
x_future=zeros(2,M,steps+1);
for i=1:M
    [~, x_temp]=ode45(@(t,x) duffing(t,x),0:time_step:steps*time_step,x(:,i));
    x_future(:,i,:)=x_temp.';
end

%% exact evolution of future coefs
fut_coefs=zeros(r,steps+1);
fut_coefs(:,1)=coefs;
for i=1:steps
    fut_coefs(:,i+1)=K2*fut_coefs(:,i);
end

%% exact errors
exact_errors = zeros(steps+1,1);
coefs_svd=U*S2*coefs;
fut_coefs_svd=U*S2*fut_coefs;

%% 
PX_temp=zeros(M,N);
pf=parfor_progress(steps);
pfcleanup=onCleanup(@() delete(pf));
G_a=zeros(N,N,steps);
G_b=zeros(N,N,steps);
for i=1:steps
    for j=1:sqrt(N)
        for k=1:sqrt(N)
            f1=chebpoly(j-1,[-2.5 2.5]);
            f2=chebpoly(k-1,[-2.5 2.5]);
            for l=1:M
                PX_temp(l,(j-1)*sqrt(N)+k)=f1(x_future(1,l,i+1))*f2(x_future(2,l,i+1));
            end
        end
    end
    G_a(:,:,i)=PX_temp'*PX/M;
    G_b(:,:,i)=PX_temp'*PX_temp/M;
    exact_errors(i+1)=sqrt(fut_coefs(:,i+1)'*fut_coefs(:,i+1)-2*real(fut_coefs_svd(:,i+1)'*G_a(:,:,i).'*coefs_svd)+coefs_svd'*G_b(:,:,i)*coefs_svd);
    parfor_progress(pf);
end

%% strict error bounds
proj_errors=kmd_error_bound_first_order_gelfand2(coefs,steps,G2,A2,L2,K2);
proj_errors2=kmd_error_bound_first_order(coefs,steps,G2,A2,L2,K2);

%%
normalization1=coefs'*coefs;

%% Plot relative forecast errors compared to test data
figure
semilogy(0:1:steps,exact_errors/sqrt(normalization1),'linewidth',2)
hold on
plot(0:1:steps,proj_errors2/sqrt(normalization1),'linewidth',2)  
plot(0:1:steps,proj_errors/sqrt(normalization1),'linewidth',2)  
box on
ax=gca; ax.FontSize=18; axis tight
title('Duffing oscillator','fontsize',18,'interpreter','latex')
xlabel('Time steps','interpreter','latex','fontsize',18)
ylabel('$L^2$ norm forecast error','interpreter','latex','fontsize',18)
legend({'Exact error','$E_1$','$E_2$'},'interpreter','latex','fontsize',16,'location','best')

%%%%%%%%%%%%%%%%%%%%%%%% L2 pad comparison %%%%%%%%%%%%%%%%%%%%%

%% Plot principal angles and vectors

[theta,U1,U2,J]=KoopAngles(G,A,L);

%% convert U1 to lie entirely inside V
[size1,size2]=size(U1);
V1 = U1(1:size1/2,:) + K*U1(size1/2+1:size1,:);

%% convert U2 to lie entirely inside KV
V2 = K\U2(1:size1/2,:) + U2(size1/2+1:size1,:);

%% 
exact_errors_pad=zeros(1,steps+1);

%% get r smallest angles
r=49;
U=V1(:,1:r);
G3=U'*G*U; A3=U'*A*U; L3=U'*L*U; K3=A3;

%% get coefficients of x variable
coefs_pad = (PX*U)\transpose(x);
coefs_pad = coefs_pad(:,1);

%% evolution
fut_coefs_pad=zeros(r,steps+1);
fut_coefs_pad(:,1)=coefs_pad;
for i=1:steps
    fut_coefs_pad(:,i+1)=K3*fut_coefs_pad(:,i); 
end

%% exact errors
for i=1:steps
    exact_errors_pad(i+1)=sqrt(fut_coefs_pad(:,i+1)'*G3*fut_coefs_pad(:,i+1)-2*real(fut_coefs_pad(:,i+1)'*U'*G_a(:,:,i).'*U*coefs_pad)+coefs_pad'*U'*G_b(:,:,i)*U*coefs_pad);
end


%% strict error bounds
[proj_errors_pad,K_op_pad]=kmd_error_bound_first_order_gelfand2(coefs_pad,steps,G3,A3,L3,K3);
proj_errors2_pad=kmd_error_bound_first_order(coefs_pad,steps,G3,A3,L3,K3);

%%
normalization2=coefs_pad'*coefs_pad;

%% Plot relative forecast errors compared to test data
figure
semilogy(0:1:steps,exact_errors/sqrt(normalization1),'linewidth',2)
hold on
plot(0:1:steps,min(proj_errors,proj_errors2)/sqrt(normalization1),'linewidth',2)  
plot(0:1:steps,exact_errors_pad/sqrt(normalization2),'linewidth',2)
plot(0:1:steps,min(proj_errors_pad,proj_errors2_pad)/sqrt(normalization2),'linewidth',2)  
box on
ax=gca; ax.FontSize=18; axis tight
title('Duffing oscillator','fontsize',18,'interpreter','latex')
xlabel('Time steps','interpreter','latex','fontsize',18)
ylabel('$L^2$ norm forecast error','interpreter','latex','fontsize',18)
legend({'Exact error','$E_1$','Exact error pad','$E_1$ pad'},'interpreter','latex','fontsize',16,'location','best')

%%%%%%%%%%%%%%%%%%%%%%%%%% L2 expected errors %%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% evals of kernel
Lambda = exp(-(0:1:(length(coefs)-1))/(1000)).';

%% get means and variance of norms
NN=10000;
K_temp=eye(length(G2));
K_op_avg=zeros(steps+1,1); K_op_avg(1)=1;
pf = parfor_progress(steps);
pfcleanup = onCleanup(@() delete(pf));
for i=1:steps
    K_temp=K2*K_temp;
    KGK=K_temp'*G2*K_temp;
    for j=1:NN
        v=randn(length(coefs),1);
        v=v.*(Lambda);
        K_op_avg(i+1) = K_op_avg(i+1) + sqrt((v'*KGK*v)/(v'*v));
    end
    K_op_avg(i+1)=K_op_avg(i+1)/NN;
    parfor_progress(pf);
end

%%
Q=L2-A2'*K2; Q=(Q+Q')/2;
Kg=zeros(length(coefs),steps); Kg(:,1)=coefs;
KgQKg=zeros(steps,1); KgQKg(1)=Kg(:,1)'*Q*Kg(:,1);
for i=2:steps
    Kg(:,i)=K2*Kg(:,i-1);
    KgQKg(i)=Kg(:,i)'*Q*Kg(:,i);
end
KgQKg=real(sqrt(KgQKg));
exp_error2=zeros(steps+1,1);
for i=1:steps
    for j=0:i-1
        exp_error2(i+1)=exp_error2(i+1)+(K_op_avg(2)^j*KgQKg(i-j))^2;
    end
    exp_error2(i+1)=sqrt(exp_error2(i+1));
end

%% Plot relative forecast errors compared to test data w averaging and std 
figure
semilogy(0:1:steps,exact_errors/normalization1,'linewidth',2)
hold on
plot(0:1:steps,proj_errors2/normalization1,'linewidth',2)  
plot(0:1:steps,exp_error2/normalization1,'linewidth',2)
box on
ax=gca; ax.FontSize=18; axis tight
title('Duffing oscillator','fontsize',18,'interpreter','latex')
xlabel('Time steps','interpreter','latex','fontsize',18)
ylabel('$L^2$ norm forecast error','interpreter','latex','fontsize',18)
legend({'Exact error','$E_1$','$E_2$'},'interpreter','latex','fontsize',16,'location','best')


%%
%%%%%%%%%%%% save data %%%%%%%%%%%%%%%%%%
save('big_figure_plots\\duffing_l2_error_bounds.mat','steps','exact_errors','proj_errors','proj_errors2','exact_errors_pad','proj_errors_pad','proj_errors2_pad','exp_error2','normalization1','normalization2')

%%%%%%% kernelized error bounds %%%%%%%%%%%%
%%
rng(0)

%% generate snapshot data
time_step=0.25;
steps=20;
no_traj=200;
M=no_traj*steps;
init=4*(rand(2,no_traj+1)-0.5); %4
x=zeros(2,M);
y=zeros(2,M);
for i=1:no_traj
    [~, x1]=ode45(@(t,z) duffing(t,z),0:time_step:(steps*time_step),init(:,i));
    x1=x1';
    for j=1:steps
        x(:,(i-1)*steps+j)=x1(:,j);
        y(:,(i-1)*steps+j)=x1(:,j+1);
    end
end
[~, x_test]=ode45(@(t,x) duffing(t,x),0:time_step:(steps*time_step),init(:,no_traj+1));
x_test=x_test';

%%
scale = 1/mean(std(x.'));

%% compute operator folding matrices
ker=@(x,t) kernel(x,t,scale);
[G,A,R]=generate_matrices_kernelized(x,y,ker);

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
x0=init(:,no_traj+1);
Kx0_vals=zeros(M,1);
for i=1:M
    Kx0_vals(i)=ker(x0,x(:,i));
end
coefs=(U*Sig\Kx0_vals); %coefficients of K_x0 observables in terms of data kernels
fut_coefs=zeros(r,steps+1);
fut_coefs(:,1)=coefs;
for i=1:steps
    fut_coefs(:,i+1)=K2*fut_coefs(:,i); %evolution of coefs under PF op
end
x_rkhs=x*U*S2*conj(fut_coefs);

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
delta=sqrt(ker(x0,x0)-coefs'*coefs);

%% strict error bounds
[proj_errors,K_op]=kmd_error_bound_first_order_gelfand2(coefs,steps,G2,A2,R2,K2);
[proj_errors2,K_op2]=kmd_error_bound_first_order(coefs,steps,G2,A2,R2,K2);
delta_errors=delta*K_op;

%% Plot relative forecast errors compared to test data w averaging and std 
figure
semilogy(0:1:steps,exact_errors,'linewidth',2)
hold on
plot(0:1:steps,(proj_errors+delta_errors),'linewidth',2)  
plot(0:1:steps,(min(proj_errors,proj_errors2)+delta_errors),'linewidth',2)  
box on
ax=gca; ax.FontSize=18;
title('Duffing oscillator','fontsize',18,'interpreter','latex')
xlabel('Time steps','interpreter','latex','fontsize',18)
ylabel('RKHS norm forecast error','interpreter','latex','fontsize',18)
legend({'Exact error','Error bound','Projection error','Delta error'},'interpreter','latex','fontsize',16,'location','best')

%%%%%%%%%%%%%%%%%%%%%%% kernel expected errors %%%%%%%%%%%%%%%%%%%%%%%%%%

%%
NN=10000;
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
        v=v.*(Lambda);
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

%% normalization
normalization=coefs'*coefs;

%% Plot relative forecast errors compared to test data w averaging and std 
figure
semilogy(0:1:steps,exact_errors/sqrt(normalization),'linewidth',2)
hold on
plot(0:1:steps,(proj_errors+delta_errors)/sqrt(normalization),'linewidth',2)  
plot(0:1:steps,(exp_error)/sqrt(normalization),'linewidth',2)  
box on
ax=gca; ax.FontSize=18;
title('Duffing oscillator','fontsize',18,'interpreter','latex')
xlabel('Time steps','interpreter','latex','fontsize',18)
ylabel('RKHS norm forecast error','interpreter','latex','fontsize',18)
legend({'Exact error','Error bound','Expected error'},'interpreter','latex','fontsize',16,'location','best')

%%
%%%%%%%%%%%% save data %%%%%%%%%%%%%%%%%%
save('big_figure_plots\\duffing_kernel_error_bounds.mat','steps','exact_errors','proj_errors','delta_errors','exp_error','delta_errors_avg','normalization')

%% define kernel
function ker=kernel(x,t,scale)
    r=vecnorm(x-t);
    r=scale*r;
    ker=zeros(1,size(x,2));
    ker(r>0)=(r(r>0)).^(2).*besselk(-2,r(r>0));
    ker(r==0)=2;
end

%% define duffing oscillator
function dxdt = duffing(~,x)
    dxdt=zeros(2,1);
    dxdt(1)=x(2);
    dxdt(2)=-0.05*x(2)+x(1)-x(1).^3;
end