function [error,K_op] = kmd_error_bound_first_order(g,n,G,A,R,K_star,unitary)
    if nargin<6
        K_star=G\A;
    end
    if nargin<7
        unitary=0;
    end
    Kg=zeros(length(g),n);
    Kg(:,1)=g;
    for i=2:n
        Kg(:,i)=K_star*Kg(:,i-1);
    end
    Q=R-A'*K_star;
    if unitary==0
        SQ=chol(G);
        K_op=norm(SQ*K_star/SQ);
    else
        K_op=1;
    end
    error=zeros(n+1,1);
    for i=1:n
        for j=0:i-1
            error(i+1)=error(i+1)+K_op^(2*j)*Kg(:,i-j)'*Q*Kg(:,i-j);
        end
    end
    error=sqrt(error);
    K_op=K_op.^(0:1:n);
    K_op=K_op.';
end