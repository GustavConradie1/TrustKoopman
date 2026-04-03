function [theta,U1,V1,J] = KoopAngles(G,A,L,B,ec)


% Default parameters
if nargin < 4
    B = eye(size(G));
    ec = 10^(-30); 
elseif nargin < 5
    ec = 10^(-30); 
end

% Fail safe for G and L Hermitian
G = (G+G')/2; L = (L+L')/2;

% Compress matrices to subspace
Gv = B'*G*B;
Av = B'*A*B;
Lv = B'*L*B;

% Form weight matrix and eval decomp
J = [Gv, Av; Av', Lv]; J = (J+J')/2;
[U,E] = eig(J,'vector'); E = E(:);

% Form range matrices
I = find(E>ec);
C1 = sqrt(E(I)).*(U(1:size(Gv,1),I))';
C2 = sqrt(E(I)).*(U(1+size(Gv,1):end,I))';
Q1 = orth(C1);
Q2 = orth(C2);

% Compute principal angles and vectors
[theta,U1,V1] = subspacea(Q1,Q2);

% Convert back to original basis
U1 = U(:,I)*(1./sqrt(E(I)).*U1);
V1 = U(:,I)*(1./sqrt(E(I)).*V1);

end

