function [error,K_op] = kmd_error_bound_first_order_gelfand2(g,n,G,A,R,K_star,unitary)
    if nargin<6
        K_star=G\A;
    end
    if nargin<7
        unitary=0;
    end
    Q=R-A'*K_star; Q=(Q+Q')/2;
    SQ=chol(G);
    K_op=ones(n+1,1);
    if unitary==0
        K_temp=eye(length(G));
        for i=1:n
            K_temp=K_star*K_temp;
            K_op(i+1)=norm(SQ*K_temp/SQ);
        end
    end
    Kg=zeros(length(g),n);
    KgQKg=zeros(n,1);
    Kg(:,1)=g;
    KgQKg(1)=Kg(:,1)'*Q*Kg(:,1);
    for i=2:n
        Kg(:,i)=K_star*Kg(:,i-1);
        KgQKg(i)=Kg(:,i)'*Q*Kg(:,i);
    end
    KgQKg=real(sqrt(KgQKg));
    error=zeros(n+1,1);
    for i=1:n
        for j=0:i-1
            error(i+1)=error(i+1)+K_op(j+1)*KgQKg(i-j);
        end
    end
end