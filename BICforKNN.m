function [bic] = BICforKNN(X,Beta,Quantile,Y)
q = size(Y,2); p = size(X,2); n = size(Y,1);
bic=0;

parfor e=1:q
set = (e-1)*p+1:e*p
[theta_e] = Beta(:,set);
Z = Y(:,e) - sum(X.*theta_e,2);
qt = (Z>0).*Quantile' + (Z<0).*(Quantile'-1);
values_lik = log(sum(Z.*qt,"all"))
theta_round = round(theta_e,2);
sparsity=[];
for ii=1:p
    sparsity=[sparsity,length(unique(theta_round(:,ii)))];
end
bic = bic + values_lik + sum(sparsity)*log(n)/n
end

end