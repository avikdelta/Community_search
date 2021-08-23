function [error permidx] = ClusterErr(Theta, Theta_Hat, kIdx)

[n, K] = size(Theta);
[n1, K1] = size(Theta_Hat);

if kIdx == 0
    if K1 < K
        Theta_Hat = double(rand(n,K)>.5);
    end

    P = perms(1:K);

    error = n*K + 10;

    for i = 1:factorial(K)
        E = Theta_Hat(:,P(i,:));
        temp = sum(sum(double(Theta ~= E),1));

        if temp < error
            error = temp;
            permidx = i;
        end

    end
else
    eArr = zeros(1,K1);
    for k = 1:K1
        eArr(k) = sum(double(Theta(:,kIdx)~=Theta_Hat(:,k)));
    end
    error = min(eArr);
    permidx = 0;
end

end