%% Generating Data
lambdaset1 = [0.5, 0.7, 0.9, 1.1, 1.3, 1.5, 1.7, 1.9];
lambdaset2 = [0.5,1];
[l1, l2] = ndgrid(lambdaset1, lambdaset2);
combinations = [l1(:), l2(:)];
n= 2000; q=5; p=6;
T = rand(1,n);
Quantile = rand(1,n);
Quantile = 0.9*Quantile + 0.05; 
errs = normrnd(0,1,[n,q]);
X = rand(2000,p);
B1 = reshape(ones(p*q,1),q,p);
g1 = sqrt(6/7)*(T+Quantile)'; 
Beta =  g1*B1(:)';
Y = zeros(n,q);
parfor i=1:n
    xi = X(i,:); Iq = eye(q);
    y = kron(xi,Iq)*(Beta(i,:))'
    Y(i,:) = y
end
Y=Y+errs;

%% KNN method
lamknn = [0.3,0.4,0.5,0.8,0.9,1.0];
thetaset2 = cell(1,6);
knn_bic = [];
parfor lam=1:6
    lambda = lamknn(lam)
    [theta1, knn_bic1] = MVC(X,Quantile,T,Y,5,lambda,500,0.1,0.5);
    thetaset2{lam} = theta1
    knn_bic = [knn_bic,knn_bic1]
end
min_knn = min(knn_bic);
min_idx = find(knn_bic==min_knn);
theta = thetaset2{min_idx};

%% RVC
BIC = [];
D_set2 = cell(1,8);
parfor tp=1:8
    tunparset = combinations(tp,:)
    lambda1 = tunparset(1)
    lambda2 = tunparset(2)
    [D,f] = RVC(X,Quantile,T,Y,theta,lambda1,lambda2,5)
    D_set2{tp} = D
    [bicvalue] = BICforRVC(X,Quantile,Y,D,f)
    BIC = [BIC, bicvalue]
end
min_BIC = min(BIC);
min_ind = find(BIC==min_BIC);
K = rank(D_set2{min_ind});
D_set = cell(1,6);
f_set = cell(1,6);
BIC = [];
parfor lam=1:6
    lambda2 = lamknn(lam)
    [D,f] = RVC(X,Quantile,T,Y,theta,0,lambda2,K);
    D_set{lam} = D
    f_set{lam} = f
    Result=f*D'
    re=1:size(X,2)*size(Y,2)
    reconst=reshape(re,q,p)'
    Result=Result(:,reconst(:))
    bic = BICforKNN(X,Result,Quantile,Y)
    BIC = [BIC, bic]
end
minBIC = find(BIC==min(BIC));
D=D_set{minBIC};
f=f_set{minBIC};

%% Reconstructing Results
re = 1:p*q;
reconst = reshape(re,q,p)';
Result = f*D';
Result = Result(:,reconst(:));

%% Visualizing Results 
color_limits=[];
for ii = 1:30
subplot(5,6,ii)
scatter(T,Quantile, 14, Result(:,ii)) 
ttt=xlabel('T'); ttt.FontSize=6;
ylabel('\tau')
title(sprintf('\\beta_%d',ii),'FontSize', 8)
uuu=colorbar;
color_limits=[color_limits;uuu.Limits]; uuu.TickLength=0.05;
end
%%

exportgraphics(gcf, 'fig_demo_simul.jpg', 'Resolution', 300);

