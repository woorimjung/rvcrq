function [bicvalue] = BICforRVC(X,Quantile,Y,D,f)

n = size(Y,1); q = size(Y,2); p = size(X,2)-1; 
K= size(f,2);
beta_bic = f*D';
    parfor i=1:n
        Xi = X(i,:)';
        bi = reshape(beta_bic(i,:),q,p+1)
        Xb(i,:) = bi*Xi
    end
Z = Y - Xb;
Q = repmat(Quantile',1,q);
qt = (Z>0).*Q + (Z<0).*(Q-1);
ch_val = log(sum(Z.*qt,"all"));

f2_round = round(f,2);
sparsity=[];
for k=1:K
    sparsity=[sparsity,length(unique(f2_round(:,k)))];
end

sval = rank(D);
bicvalue = ch_val + (log(n)/n)*sum(sparsity) + (log(n)/n)*sval*(p+1)*q;
end