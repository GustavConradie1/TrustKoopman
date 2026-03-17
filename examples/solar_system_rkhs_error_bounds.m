%% Koopman/PF analysis of antarctic sea ice concentration data
clear
rng(0)

%% load data sets
load('states_pluto.mat')
load('states_charon.mat')
load('states_hydra.mat')
load('states_styx.mat')
load('states_nix.mat')
load('states_kerberos.mat')

DATA=[states_pluto states_charon states_hydra states_styx states_nix states_kerberos];
DATA = DATA(1:365*80,:);
DATA=DATA.';
[d,N]=size(DATA);
no_bodies=d/6;
delta_t=10;

%%
means = zeros(d,1);
stds = zeros(d,1);
for i=1:d
    means(i)=mean(DATA(i,: ));
    stds(i)=std(DATA(i,:));
    DATA(i,:)=(DATA(i,:)-means(i))/stds(i);
end

%% delay embedding to get acceleration esque data 
steps=72;
maxx=(N-delta_t*steps);
minn=1;
x=DATA(:,minn:maxx-delta_t);
y=DATA(:,minn+delta_t:maxx);
M=maxx-minn-delta_t+1;

%% compute matrices
s_vals=1/10000*ones(no_bodies,1);
v_vals=1*ones(no_bodies,1);
ker=@(x,t) kernel_matern_product(x,t,s_vals,v_vals,no_bodies);
[G,A,R]=generate_matrices_kernelized(x,y,ker);

%% change of basis
r=5000; 
tic
[U,Sig]=eigs(G,r,'largestabs');
toc
Sig = sqrt(Sig); S2 = diag(1./diag(Sig));
A2 = S2*(U'*A*U)*S2; 
R2 = S2*(U'*R*U)*S2; R2 = (R2+R2')/2;
G2 = eye(r);
K2 = A2;

%% predict future trajectory
x0=DATA(:,maxx);
Kx0_vals=zeros(M,1);
for i=1:M
    Kx0_vals(i)=ker(x0,x(:,i));
end
coefs=(U*Sig\Kx0_vals); 
fut_coefs=zeros(r,steps+1);
fut_coefs(:,1)=coefs;
for i=1:steps
    fut_coefs(:,i+1)=K2*fut_coefs(:,i);
end
x_rkhs=x*U*S2*conj(fut_coefs);

%%
real_data=DATA(:,maxx+(0:delta_t:delta_t*steps));
exact_errors=sqrt(sum(abs(x_rkhs(end-d+1:end,:)-real_data(end-d+1:end,:)).^2,1)./sum(abs(real_data(end-d+1:end,:)).^2,1));

%% Plot relative forecast errors
figure
plot(0:delta_t:delta_t*steps,exact_errors,'linewidth',2)

grid on
legend('Exact','interpreter','latex','fontsize',18,'location','best')
xlabel('Time steps (days)','interpreter','latex','fontsize',18)
ylabel('Relative pointwise errors','interpreter','latex','fontsize',18)
title('Pluto-Charon system','interpreter','latex','fontsize',18)

%% plot exact versus predicted trajectories
for i=19:24
    figure
    box on
    hold on
    if mod(i,6)==1 || mod(i,6)==2 || mod(i,6)==3
        p1=plot(0:delta_t:delta_t*steps,stds(i)*x_rkhs(i,:)+means(i)-(stds(mod(i,6))*x_rkhs(mod(i,6),:)+means(mod(i,6))),'linewidth',2,'color',[0.9290    0.6940    0.1250]);
        p2=plot(0:delta_t:delta_t*steps,stds(i)*real_data(i,:)+means(i)-(stds(mod(i,6))*real_data(mod(i,6),:)+means(mod(i,6))),'--','linewidth',2,'color',[0    0.4470    0.7410]);
    else
        p1=plot(0:delta_t:delta_t*steps,stds(i)*x_rkhs(i,:)+means(i),'linewidth',2,'color',[0.9290    0.6940    0.1250]);
        p2=plot(0:delta_t:delta_t*steps,stds(i)*real_data(i,:)+means(i),'--','linewidth',2,'color',[0    0.4470    0.7410]);
    end
    legend([p2 p1],'Exact','kEDMD','interpreter','latex','fontsize',18,'location','best')
    xlabel('Time steps (days)','interpreter','latex','fontsize',18)
    ylabel('Variable','interpreter','latex','fontsize',18)
    ax=gca; ax.FontSize=18; xlim('tight')
    if mod(i,6)==1
        title('Styx, $x$-coordinate','interpreter','latex','fontsize',18)
    elseif mod(i,6)==2
        title('Styx, $y$-coordinate','interpreter','latex','fontsize',18)
    elseif mod(i,6)==3
        title('Styx, $z$-coordinate','interpreter','latex','fontsize',18)
    elseif mod(i,6)==4
        title('Styx, $\dot{x}$-coordinate','interpreter','latex','fontsize',18)
    elseif mod(i,6)==5
        title('Styx, $\dot{y}$-coordinate','interpreter','latex','fontsize',18)
    elseif mod(i,6)==0
        title('Styx, $\dot{z}$-coordinate','interpreter','latex','fontsize',18)
    end
end

%%
G_extra=zeros(steps+1,M);
for i=1:steps+1
    for j=1:M
        G_extra(i,j)=ker(real_data(:,i),x(:,j));
    end
end

%%
exact_errors_rkhs=zeros(steps+1,1);
for i=1:steps+1
    exact_errors_rkhs(i)=sqrt(ker(real_data(:,i),real_data(:,i))+fut_coefs(:,i)'*fut_coefs(:,i)-2*G_extra(i,:)*U*S2*fut_coefs(:,i));
end

%% error bounds on evolution of K_{x_0}
delta=sqrt(ker(x0,x0)-coefs'*coefs);

%% strict error bounds
tic
[proj_errors,K_op]=kmd_error_bound_first_order_gelfand2(coefs,steps,G2,A2,R2,K2);
[proj_errors2,K_op2]=kmd_error_bound_first_order(coefs,steps,G2,A2,R2,K2);
toc
delta_errors=delta*K_op;

%%
NN=100;
Lambda = exp(-(0:1:(length(coefs)-1))/(1000)).';

%% get means and variance of norms
K_temp=eye(r);
K_op_avg=zeros(steps+1,1); K_op_avg(1)=1;
K_op_avg2=zeros(steps+1,1); K_op_avg2(1)=1;
K_op_std=zeros(steps+1,1);
pf = parfor_progress(steps);
pfcleanup = onCleanup(@() delete(pf));
for i=1:steps
    K_temp=K2*K_temp;
    KGK=K_temp'*G2*K_temp;
    for j=1:NN
        v=randn(r,1);
        v=v.*Lambda;
        K_op_avg(i+1) = K_op_avg(i+1) + sqrt((v'*KGK*v)/(v'*v));
        K_op_avg2(i+1) = K_op_avg2(i+1) + (v'*KGK*v)/(v'*v);
    end
    K_op_avg(i+1)=K_op_avg(i+1)/NN;
    K_op_avg2(i+1)=K_op_avg2(i+1)/NN;
    K_op_std(i+1)=sqrt(K_op_avg2(i+1)-K_op_avg(i+1)^2);
    parfor_progress(pf);
end

%%
k=3;
delta_errors_avg=delta*K_op_avg;
delta_errors_above=delta*(K_op_avg+k*K_op_std);
delta_errors_below=delta*(K_op_avg-k*K_op_std);

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
exp_error2=exp_error2+delta_errors_avg;
exp_error2_above=exp_error2_above+delta_errors_above;
exp_error2_below=exp_error2_below+delta_errors_below;

%%
exp_error=zeros(steps+1,1); exp_error(1)=delta;
exp_error_above=zeros(steps+1,1); exp_error_above(1)=delta;
exp_error_above2=zeros(steps+1,1); exp_error_above2(1)=delta^2;
exp_error_below=zeros(steps+1,1); exp_error_below(1)=delta;
exp_error_below2=zeros(steps+1,1); exp_error_below2(1)=delta^2;
pf = parfor_progress(steps);
pfcleanup = onCleanup(@() delete(pf));
for i=1:steps
    temp=zeros(1,i+1);
    temp_above=zeros(1,i+1);
    temp_below=zeros(1,i+1);
    temp(i+1)=delta_errors_avg(i+1);
    temp_above(i+1)=delta_errors_above(i+1);
    temp_below(i+1)=delta_errors_below(i+1);
    for j=0:i-1
        temp(j+1)=K_op_avg(j+1)*KgQKg(i-j);
        temp_above(j+1)=(K_op_avg(j+1)+k*K_op_std(j+1))*KgQKg(i-j);
        temp_below(j+1)=(K_op_avg(j+1)-k*K_op_std(j+1))*KgQKg(i-j);
    end
    for j=1:NN
        v=randn(length(coefs),i+1);
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
    exp_error(i+1)=exp_error(i+1)/NN; exp_error_above(i+1)=exp_error_above(i+1)/NN; exp_error_above2(i+1)=exp_error_above2(i+1)/NN;
    exp_error_below(i+1)=exp_error_below(i+1)/NN; exp_error_below2(i+1)=exp_error_below2(i+1)/NN;
    parfor_progress(pf);
end
exp_error_above_std=sqrt(exp_error_above2-exp_error_above.^2); exp_error_below_std=sqrt(exp_error_below2-exp_error_below.^2);
exp_error_above=exp_error_above+k*exp_error_above_std; exp_error_below=exp_error_below-k*exp_error_below_std;

%%
normalization=sqrt(coefs'*coefs);
colors=orderedcolors('gem');
fontsize=36;
axissize=28;
linesize=5;

%% Plot relative forecast errors compared to test data w averaging and std 
figure
semilogy(0:delta_t:delta_t*steps,exact_errors_rkhs/normalization,'linewidth',linesize,'color','black')
hold on
plot(0:delta_t:delta_t*steps,(proj_errors+delta_errors)/normalization,'linewidth',linesize,'color',1.2*colors(1,:))  
plot(0:delta_t:delta_t*steps,exp_error/normalization,'linewidth',linesize,'color',colors(2,:))  
patch([((0:delta_t:delta_t*steps).');flip((0:delta_t:delta_t*steps).')],real([(exp_error_below);flip(exp_error_above)])/normalization,1*colors(2,:),'FaceAlpha',.2)
box on
ax=gca; ax.FontSize=axissize; axis tight
title('Pluto-Charon system','fontsize',fontsize,'interpreter','latex')
xlabel('Time steps (days)','interpreter','latex','fontsize',fontsize)
ylabel('RKHS Error','interpreter','latex','fontsize',fontsize)
legend({'Exact','Bound','Expected$\;$'},'interpreter','latex','fontsize',fontsize,'location','best')

%%
%%%%%%%%%%%% save data %%%%%%%%%%%%%%%%%%
save('big_figure_plots\\solar_system_var_results.mat','delta_t','steps','proj_errors','delta_errors','exp_error','exp_error_below','exp_error_above','real_data','exact_errors_rkhs','exact_errors','x_rkhs','means','stds','normalization')
disp('saved')


%%
function ker=kernel_matern_product(x,t,s,v,no)
    ker=ones(1,size(x,2));
    eps=10^(-16);
    for j=1:no
        r=vecnorm(x((j-1)*6+1:j*6,:)-t((j-1)*6+1:j*6,:));
        ker(r>0)=ker(r>0).*(s(j)*r(r>0)).^(v(j)).*besselk(v(j),s(j)*r(r>0));
        ker(r==0)=ker(r==0).*(s(j)*eps).^(v(j)).*besselk(v(j),s(j)*eps);
    end
end