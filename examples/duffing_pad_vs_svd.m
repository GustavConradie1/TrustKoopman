clear
% needs chebfun to run

%% generate snapshot data
rng(0)
time_step=0.25;
steps=20;
test_steps=20;
no_traj=100;
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

%% get coefficients of full state observable for full model
coefs = PX\transpose(x);
coefs=coefs(:,1); 
normalization=sqrt(coefs'*G*coefs);

%% evolution
fut_coefs=zeros(N,test_steps+1);
fut_coefs(:,1)=coefs;
for i=1:test_steps
    fut_coefs(:,i+1)=K*fut_coefs(:,i);
end

%% exact errors of observables
x_future=zeros(2,M,steps+1);
for i=1:M
    [~, x_temp]=ode45(@(t,x) duffing(t,x),0:time_step:steps*time_step,x(:,i));
    x_future(:,i,:)=x_temp.';
end

%% generate G_a and G_b
G_a=zeros(N,N,test_steps);
G_b=zeros(N,N,test_steps);

PX_temp=zeros(M,N);
pf=parfor_progress(test_steps);
pfcleanup=onCleanup(@() delete(pf));
for i=1:test_steps
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
    parfor_progress(pf);
end

%%
exact_errors=zeros(test_steps+1,1);
for i=1:test_steps
    exact_errors(i+1)=sqrt(fut_coefs(:,i+1)'*G*fut_coefs(:,i+1)-2*real(fut_coefs(:,i+1)'*G_a(:,:,i).'*coefs)+coefs'*G_b(:,:,i)*coefs);
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
idx=find(theta<0.01); r=length(idx);
U=V1(:,idx);
G2=U'*G*U; A2=U'*A*U; L2=U'*L*U; K2=A2;

%% coefs 
coefs = (PX*U)\transpose(x);
coefs = coefs(:,1);
normalization_pad=sqrt(coefs'*coefs);

%% evolution
fut_coefs=zeros(r,test_steps+1);
fut_coefs(:,1)=coefs;
for i=1:test_steps
    fut_coefs(:,i+1)=K2*fut_coefs(:,i);
end

%%
exact_errors_pad=zeros(test_steps+1,1);
for i=1:test_steps
    exact_errors_pad(i+1)=sqrt(fut_coefs(:,i+1)'*G2*fut_coefs(:,i+1)-2*real(fut_coefs(:,i+1)'*U'*G_a(:,:,i).'*U*coefs)+coefs'*U'*G_b(:,:,i)*U*coefs);
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
coefs = (PX*U*S2)\transpose(x);
coefs = coefs(:,1);
normalization_svd=sqrt(coefs'*coefs);

%% evolution
fut_coefs=zeros(r,test_steps+1);
fut_coefs(:,1)=coefs;
for i=1:test_steps
    fut_coefs(:,i+1)=K3*fut_coefs(:,i); %evolution of coefs under PF op (transpose?)
end

%%
exact_errors_svd=zeros(test_steps+1,1);
for i=1:test_steps
    exact_errors_svd(i+1)=sqrt(fut_coefs(:,i+1)'*G3*fut_coefs(:,i+1)-2*real(fut_coefs(:,i+1)'*S2*U'*G_a(:,:,i).'*U*S2*coefs)+coefs'*S2*U'*G_b(:,:,i)*U*S2*coefs);
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
r_range=5:1:N;
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
    coefs = (PX*U)\transpose(x);
    coefs = coefs(:,1);
    normalization_pad=sqrt(coefs'*coefs);

    %% evolution
    fut_coefs=zeros(r,test_steps+1);
    fut_coefs(:,1)=coefs;
    for i=1:test_steps
        fut_coefs(:,i+1)=K2*fut_coefs(:,i); 
    end

    %% exact errors
    for i=1:test_steps
        error_pad(ii,i)=sqrt(fut_coefs(:,i+1)'*G2*fut_coefs(:,i+1)-2*real(fut_coefs(:,i+1)'*U'*G_a(:,:,i).'*U*coefs)+coefs'*U'*G_b(:,:,i)*U*coefs)/normalization_pad;
    end

    %% truncate via usual SVD
    [U,Sig]=eigs(G,r,'largestabs');
    Sig = sqrt(Sig); S2 = diag(1./diag(Sig));
    A3 = S2*(U'*A*U)*S2; 
    L3 = S2*(U'*L*U)*S2; L3 = (L3+L3')/2;
    G3 = eye(r);
    K3 = A3;
    
    %% coefs 
    coefs = (PX*U*S2)\transpose(x);
    coefs = coefs(:,1);
    normalization_svd=sqrt(coefs'*coefs);
    
    %% evolution
    fut_coefs=zeros(r,test_steps+1);
    fut_coefs(:,1)=coefs;
    for i=1:test_steps
        fut_coefs(:,i+1)=K3*fut_coefs(:,i); 
    end
        
    %%
    for i=1:test_steps
        error_svd(ii,i)=sqrt(fut_coefs(:,i+1)'*G3*fut_coefs(:,i+1)-2*real(fut_coefs(:,i+1)'*S2*U'*G_a(:,:,i).'*U*S2*coefs)+coefs'*S2*U'*G_b(:,:,i)*U*S2*coefs)/normalization_svd;
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

%%
%%%%%%%%%%%% save data %%%%%%%%%%%%%%%%%%
save('big_figure_plots\\duffing_pad_vs_svd.mat','test_steps','exact_errors','exact_errors_svd','exact_errors_pad','r_range','error_pad','error_svd')

%% define duffing oscillator
function dxdt = duffing(~,x)
    dxdt=zeros(2,1);
    dxdt(1)=x(2);
    dxdt(2)=-0.05*x(2)+x(1)-x(1).^3;
end

