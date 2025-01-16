function [theta,bic] = MVC(X,Quantile,T,Y,Neighbour,lambda,Niter, tol, eta)

q = size(Y,2); p = size(X,2); n = size(Y,1);
theta=[]; bic=0;
parfor e=1:q
[theta_e] = VC_qt_knn_admm(X,Quantile,T,Y(:,e),Neighbour,lambda,Niter, tol, eta);
Z = Y(:,e) - sum(X.*theta_e,2);
qt = (Z>0).*Quantile' + (Z<0).*(Quantile'-1);
values_lik = log(sum(Z.*qt,"all"))
theta_round = round(theta_e,2);
sparsity=[];
for ii=1:p
    sparsity=[sparsity,length(unique(theta_round(:,ii)))];
end
theta = [theta, theta_e]
bic = bic + values_lik + sum(sparsity)*log(n)/n
end

end