clear
%needs chebfun to run
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

%% Plot principal angles and vectors

[theta,U1,U2,J]=KoopAngles(G,A,L);
figure
scatter(1:N,theta,'filled')
title('$\theta$','interpreter','latex','fontsize',18)

%% convert U1 to lie entirely inside V
[size1,size2]=size(U1);
V1 = U1(1:size1/2,:) + K*U1(size1/2+1:size1,:);

%% convert U2 to lie entirely inside KV
V2 = K\U2(1:size1/2,:) + U2(size1/2+1:size1,:);

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
initial_errors=zeros(size2,1);

%%
for ii=1:size2
    %% coefs
    coefs = V1(:,ii);
    fut_coefs=K*coefs;
    
    initial_errors(ii)=sqrt(fut_coefs'*G*fut_coefs-2*real(fut_coefs'*G_a(:,:,1).'*coefs)+coefs'*G_b(:,:,1)*coefs);
end

%% best piecewise linear 2 components fits of lines
[coef1,bkpt1,~]=fitBogartz(log(theta(2:end)),log(initial_errors(2:end)),0); 
r1=linspace(theta(2),exp(bkpt1(1)),100); 
s1=linspace(exp(bkpt1(1)),theta(end),100);

%% plot
colors=orderedcolors('gem');
figure
hold on
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
scatter(theta(2:end),initial_errors(2:end),15,'filled')
plot(r1,exp(coef1(1)*log(r1)+coef1(2)),'linewidth',3,'color',colors(1,:))
plot(s1,exp(coef1(3)*log(s1)+coef1(4)),'linewidth',3,'color',colors(1,:))
box on
ax=gca; ax.FontSize=18; axis tight
xlabel('Principal angle','interpreter','latex','fontsize',18)
ylabel('One-step error in principal observable','interpreter','latex','fontsize',18)
exportgraphics(gcf,'duffing_results\\duffing_angles_vs_errors.pdf','ContentType','vector','BackgroundColor','none')

%%
%%%%%%%%%%%% save data %%%%%%%%%%%%%%%%%%
save('big_figure_plots\\duffing_principal_angle.mat','theta','initial_errors')

%% define duffing oscillator
function dxdt = duffing(~,x)
    dxdt=zeros(2,1);
    dxdt(1)=x(2);
    dxdt(2)=-0.05*x(2)+x(1)-x(1).^3;
end

