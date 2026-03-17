clear
load('Re7500.mat')

%%
TT = 100;

N=100;      % number of basis functions used
M1=2*N;     % number of snapshots to compute the basis
M2=5*N;    % number of snapshots used for ResDMD matrices
M3=5*N+TT;    % number of snapshots used for ResDMD matrices

use_DMD=1;

ind1=1:M1;
ind2=(1:M3);

%% Apply ResDMD
if use_DMD~=1
    [PX,PY] = kernel_ResDMD(DATA(:,ind1),DATA(:,ind1+1),DATA(:,ind2),DATA(:,ind2+1),'N',N,'Parallel','off','cut_off',0);
else
    [~,S,V]=svd(transpose(DATA(:,ind1))/sqrt(M1),'econ');
    PX=transpose(DATA(:,ind2))*V(:,1:N)*diag(1./(diag(S(1:N,1:N))));
    PY=transpose(DATA(:,ind2+1))*V(:,1:N)*diag(1./(diag(S(1:N,1:N))));
    clear S V
end

%%
G = (PX(1:M2,:)'*PX(1:M2,:))/M2;
A = (PX(1:M2,:)'*PY(1:M2,:))/M2;
L = (PY(1:M2,:)'*PY(1:M2,:))/M2;
K=G\A;

%%%%%%%%% principal angles and observables %%%%%%%%%%%%%%%%%5

%% Principal angle and vector computation
[theta,U1,U2] = KoopAngles(G,A,L);
[size1,size2]=size(U1);
V1 = U1(1:size1/2,:) + K*U1(size1/2+1:size1,:);

%% get principal modes
modes = ([PX(1:M2,:)]*V1)\transpose(DATA(:,1:M2));

%% plot angles
figure
semilogy(theta,'.','markersize',20)
ax = gca; ax.FontSize = 16;
xlabel('$j$','interpreter','latex','fontsize',16)
title('$\theta_j(\mathcal{V},\mathcal{K}\mathcal{V})$','interpreter','latex','fontsize',20)

%% plot modes
idx=5;
u=modes(idx,:);
w = u.*abs(u).^(-0.5);
[~,I2]=sort(real(w),'ascend');
figure
scatter(xpt(I2),ypt(I2),3+2*(xpt(I2)>1),real(w(I2)),'filled');
set(gca,'xticklabel',{[]})
set(gca,'yticklabel',{[]})
colormap(coolwarm); 
s = std( real(w( (xpt(:)<1.5)&(abs(u(:))>0) ) ));
m = mean(real(w( (xpt(:)<1.5)&(abs(u(:))>0)  )));
clim([m-3*s,m+3*s])
ax=gca; ax.FontSize=18; axis equal;
grid off
axis off

%%
snapshot = DATA(:,5000).';
save('big_figure_plots\\cavity.mat','theta','modes','xpt','ypt','snapshot')

