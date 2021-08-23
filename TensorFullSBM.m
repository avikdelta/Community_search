function T = TensorFullSBM(X1,X2,X3,lambdaArr,phiMat)

[d1 N1] = size(X1);
[d2 N2] = size(X2);
[d3 N3] = size(X3);

N = min([N1,N2,N3]);

T = zeros(d1,d2,d3);
for s = 1:N
    Temp = zeros(d1,d2,d3);
    for i = 1:d3
        Temp(:,:,i) = X3(i,s)*X1(:,s)*X2(:,s)';
    end
    T = T + Temp;
end
T = T/N;

% Deflation
M = length(lambdaArr);
if M>0
    T2 = zeros(d1,d2,d3);
    for m = 1:M
        Temp = zeros(d1,d2,d3);
        for i = 1:d3
            Temp(:,:,i) = phiMat(i,m)*phiMat(:,m)*phiMat(:,m)';
        end
        T2 = T2 + lambdaArr(m)*Temp;
    end
    T = T - T2;
end


end