function[D] = upd_D(X,Y,f,Quantile,lambda1,eta_tilde,D_ini,K) 
n=size(X,1); p=size(X,2)-1; q=size(Y,2); %K =(p+1)*q; 
%K=5;
    % initial D, D_tilde, W
    D = D_ini;
    D_tilde = D;
    W = normrnd(0,1,size(D,1),K);

    % initial Z, U_tilde
    beta_upD = f*D';
    %Z = zeros(n,q);
    parfor i=1:n
        Xi = X(i,:)';
        bi = reshape(beta_upD(i,:),q,p+1)
        Xb(i,:) = bi*Xi
    end
    Z = Y - Xb;
    U_tilde = zeros(n,q);


    itrD = 0;
    errorD = 100;
    while errorD > 0.01 & itrD <= 500
        itrD = itrD+1;
        D_old = D;

         % update D
        [U, S, V] = svd(D_tilde + W/eta_tilde);
        singularvalues = diag(S); 
        singularvalues = singularvalues - lambda1/eta_tilde;
        singularvalues = max(singularvalues,0);
        S = diag(singularvalues);
        set = 1:length(singularvalues);
        D = U(:,set)*S*V';

        % Update D_tilde
        Iq = speye(q);
        UUT = zeros(K*(p+1), K*(p+1));
        small_XYUZ = zeros(K*q*(p+1),1);
        parfor i=1:n
            fi = f(i,:); Xi = X(i,:); Yi = Y(i,:); U_tilde_i = U_tilde(i,:); Zi = Z(i,:);
            fi_Xi = kron(fi,Xi);
            UUT = UUT + fi_Xi'*fi_Xi;
            v_i = (Yi-U_tilde_i/eta_tilde-Zi)';
            u_i = v_i*fi_Xi;
            small_XYUZ = small_XYUZ + u_i(:);
        end
        A_small = eye(K*(p+1)) + UUT;           
        rhs = D(:) - W(:)/eta_tilde + small_XYUZ;
        D_tilde(:) = kron(A_small \ eye(K*(p+1)), Iq) * rhs;


        % Update Z
        beta_D_tilde = f*D_tilde';
        parfor i=1:n
            Xi = X(i,:)';
            bi = reshape(beta_D_tilde(i,:),q,p+1)
            Xb(i,:) = bi*Xi
        end
        Z = Y - Xb - U_tilde/eta_tilde;
        zz = zeros(n,q);
        reQ = repmat(Quantile', 1, q);
        cond1 = n*eta_tilde*Z > reQ;
        Z(cond1) = Z(cond1) - reQ(cond1)/(n*eta_tilde);
        cond2 = n*eta_tilde*Z < reQ-1 & ~cond1;
        Z(cond2) = Z(cond2) +(1-reQ(cond2))/(n*eta_tilde);
        cond3 = n*eta_tilde*Z < reQ & n*eta_tilde*Z > reQ-1 & ~cond1 & ~cond2;
        if cond3 ~= zeros(n,q)
        Z(cond3) = zz(cond3);
        end

        U_tilde = U_tilde - eta_tilde*(Y-Xb-Z);
        % update W
        W = W - eta_tilde*(D - D_tilde);
        
        %stopping criterion
        errorD = max(norm(D-D_tilde,"fro"),eta_tilde*norm(D-D_old,"fro"));
        [errorD]
    end

end

