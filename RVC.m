function [D,f] = RVC(X,Quantile,T,Y,theta,lambda1,lambda2,K)
Neighbour=5;tolerance=0.001;Niter=100;

n=size(X,1); p=size(X,2)-1; q=size(Y,2); %K = (p+1)*q; 
%K=5;

Iq = eye(q); 

re = reshape(1:(p+1)*q,p+1,q)';
reconst = re(:);
theta1 = theta(:,reconst);
[U, S, V] = svd(theta1);
set = 1:K;
U1 = U(:,set);
S1 = S(set,:);
f = U1.*sqrt(2000); D = S1*V'./sqrt(2000); D=D';


D_ini = D;
g = f;
u = normrnd(0,0.01,n,K);

eta_tilde = 0.5; eta = 0.5;
epsilon = 3; 

error = 100;
itr = 0;
%% Reduced Varying Coefficient of Regional Quantile for Multiple Responses

while error > tolerance & itr <= Niter 
    itr = itr+1;
    D0 = D;f0 = f;
  
    % update D
    D = upd_D(X,Y,f,Quantile,lambda1,eta_tilde, D0,K);
    
    
    a = f./sqrt(n);
    % update f
    grad_a = zeros(n,K); 
    b = f*D';
        parfor i=1:n
            Xi = X(i,:)';
            bi = reshape(b(i,:),q,p+1)
            Xb(i,:) = bi*Xi
        end
    A = Y - Xb;
    Q = repmat(Quantile',1,q);
    qt = (A>0).*Q + (A<0).*(Q-1) + (A==0).*0;
    F = sum(A.*qt,"all")/n;
    parfor i=1:n
    Xi = X(i,:); 
    %Vi = kron(Xi,Iq)*D; 
    grad_a_i = zeros(1,K); ai = a(i,:); gi = g(i,:); ui = u(i,:);
    qt_i = qt(i,:);
        for k=1:K
        %Vik = Vi(:,k);
        Vik = zeros(1,q);
        Dk = D(:,k);
        for l = 1:q
            Dk_vec = Dk(l:q:end);   
            Vik(l) = Dk_vec' * Xi';
        end
        grad_a_i(k) = -(1/sqrt(n))*(qt_i*Vik') + sqrt(n)*eta*(ai(k).*sqrt(n) - gi(k) - ui(k)/eta)
        end
    grad_a(i,:) = grad_a_i;
    end

    grdt = (eye(n) - a*a')*grad_a + a*(1/2)*(a'*grad_a - grad_a'*a); %grad projection onto tangent space
    rho=1; c=(1/2)^rho; 
    [Q,R] = qr(a - c*grdt,"econ"); 
    a_t = Q;

    fnew = a_t.*sqrt(n);
    b = fnew*D';
        parfor i=1:n
            Xi = X(i,:)';
            bi = reshape(b(i,:),q,p+1)
            Xb(i,:) = bi*Xi
        end
    C = Y - Xb;
    Q = repmat(Quantile',1,q);
    qt = (C>0).*Q + (C<0).*(Q-1) + (C==0).*0;
    F_new = sum(C.*qt,"all")/n;

    parfor k=1:K 
        gk = g(:,k); uk = u(:,k);
        F = F + (eta/2) * norm(a(:,k).*sqrt(n) - gk - uk/eta)^2;
        F_new = F_new + (eta/2) * norm(a_t(:,k).*sqrt(n) - gk - uk/eta)^2;
    end
    norm_grad = norm(grdt, "fro"); 
    
    % Armijo backfitting method
    while (F_new > F - epsilon*c*norm_grad^2) & (rho<100)
        rho = rho +1;
        c = (1/2)^rho;
        [Q,R] = qr(a - c*grdt,"econ");
        a_t = Q;
    
        fnew = a_t.*sqrt(n);
        b = fnew*D';
            parfor i=1:n
                Xi = X(i,:)';
                bi = reshape(b(i,:),q,p+1)
                Xb(i,:) = bi*Xi
            end
        C = Y - Xb;
        Q = repmat(Quantile',1,q);
        qt = (C>0).*Q + (C<0).*(Q-1) + (C==0).*0;
        F_new = sum(C.*qt,"all")/n;
        for k=1:K
            F_new = F_new + (eta/2)*norm(a_t(:,k).*sqrt(n) - g(:,k) - u(:,k)/eta)^2;
        end
       [F, F_new]
    end
 
    f = a_t.*sqrt(n); 
    for kk=1:K
        if corr(f0(:,kk),f(:,kk)) < -0.5
            f(:,kk) = -f(:,kk);
        end
    end

   % update g
    Xk = [Quantile; T];
    parfor k=1:K
        yk = f(:,k) - u(:,k)/eta;
        gk = knnfl(Xk, yk', Neighbour, lambda2/eta);
        g(:,k) = gk;
    end
    

    % update u
    u = u - eta*(f - g);
    
    %stopping criterion
    error = max(max(abs(D*f' - D0*f0')));
    [error, sqrt(sum(sum((D*f'-D0*f0').^2)))]

end
end