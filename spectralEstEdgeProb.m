function [p, q] = spectralEstEdgeProb(A)

p = 0;
q = 0;
n = max(size(A));

B = A-eye(n);
S = real(eig(B));
s = sort(S,'descend');
if length(s) > 1
    s1 = s(2:end);
    diff = calcSingularDiff(s1);
    modeidx = find(diff==max(diff), 1 );
    khat = 1+modeidx;
    lambda1 = s(1);
    lambda2 = s(2);
    p = (khat*lambda1+(n-khat)*lambda2)/(n*(khat-1));
    q = (lambda1-lambda2)/n;
else
    lambda1 = s;
    lambda2 = 0;
    khat = 2;
    p = (khat*lambda1+(n-khat)*lambda2)/(n*(khat-1));
    q = (lambda1-lambda2)/n;
end

p = min(max(p,0),1);
q = min(max(q,0),p);



