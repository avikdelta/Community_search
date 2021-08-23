function [thetaMatHat, runtime] = TensorSBM(A,K,L,NumIter)
% Community search algorithm via tensor decomposition
tic;

n = size(A,1);

% Make data
X1 = A/sqrt(n);
X2 = A/sqrt(n);
X3 = A/sqrt(n);

% Compute whitening matrix
[U1 D1 V1] = svds(X1,K);
%[U2 D2 V2] = svds(X2,K);
%[U3 D3 V3] = svds(X3,K);

U2 = U1;
D2 = D1;
V2 = V1;

U3 = U1;
D3 = D1;
V3 = V1;

W1 = U1*inv(D1);
W2 = U2*inv(D2);
W3 = U3*inv(D3);

R12 = V2'*V1;
R13 = V3'*V1;

% Robust tensor decomposition
[lambdaArr, phiMat] = RobustTensorPowerSBM(X1,X2,X3,W1,W2,W3,R12,R13,K,L,NumIter);

% Recover mu
muMatHat = zeros(n,K);
for k = 1:K
    muMatHat(:,k) = lambdaArr(k)*U1*D1*phiMat(:,k);
end
alphaArrHat = 1./(lambdaArr.^2);

runtime = toc;

% Recover communities
thetaMatHat = zeros(n,K);
for i = 1:n
    %prod = A(i,:)*muMatHat;
    %kidx = find(prod==max(prod),1);
    
    row = muMatHat(i,:);
    kidx = find(row==max(row),1);
    thetaMatHat(i,kidx) = 1;
end


end