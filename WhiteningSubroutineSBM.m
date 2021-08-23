function [mu1, alpha1] = WhiteningSubroutineSBM(A1,A2,B,m1,k)

[n1 l1] = size(A1);
[n2 l2] = size(A2);

% Compute k-svd
[U1 D1 V1] = svds(A1,k);
%[U2 D2 V2] = svds(A2,k);
U2 = U1;
D2 = D1;
V2 = V1;

% Compute whitening matrix
W1 = U1*inv(D1);
W2 = U2*inv(D2);

% Whiten
Z = W1'*B*W2;

% Extract
[u1 w1 v1] = svds(Z,1);
w = U1*D1*u1;
a = u1'*W1'*m1;

mu1 = w/a;
alpha1 = a^2;

end