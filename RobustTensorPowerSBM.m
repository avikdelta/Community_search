function [lambdaArr, phiMat] = RobustTensorPowerSBM(X1,X2,X3,W1,W2,W3,R12,R13,K,L,NumIter)

% Init
[d1 N1] = size(X1);
[d2 N2] = size(X2);
[d3 N3] = size(X3);
N = min([N1,N2,N3]);

lambdaArr = [];
phiMat = [];

avgPowerIter = 0;
epsilon = 1e-3; % threshold for power iteration convergence

for k = 1:K
    disp(['Tensor power iteration for k = ' num2str(k) ', Avg. power iter. = ' num2str(avgPowerIter)]);
    
    T = TensorFullSBM(W1'*X1,R12'*W2'*X2,R13'*W3'*X3,lambdaArr,phiMat);
    %T

    % Iterate over L random starts
    lambdaPerStart = zeros(1,L);
    phiPerStart = zeros(K,L);
    powerIterCounts = zeros(1,L);
    for tau = 1:L
        % Choose theta unifromly over unit sphere
        phi = randn(K,1);
        phi = phi/norm(phi);
        % Power iterations
        diff = zeros(1,NumIter);
        for t = 1:NumIter
            % Compute power iteration update
            phiNew = zeros(K,1);        
            for i = 1:K
                temp = 0;
                for j = 1:K
                    for k = 1:K
                        temp = temp + phi(j)*phi(k)*T(i,j,k);
                    end
                end
                phiNew(i) = temp;
            end 
            % Update theta
            phiNew = phiNew/norm(phiNew);
            diff(t) = norm(phi-phiNew);
            phi = phiNew;
            
            % Check for convergence after 5 iterations
            if t>5
                % Check if maximum error over last 3 iterations is less than
                % threshold
                if max(diff(t-2:t)) < epsilon
                    powerIterCounts(tau) = t;
                    break;
                end
            end
            
        end % end of power iterations  
        % Compute lambda = M3(W*theta,W*theta,W*theta)
        temp2 = 0;
        for i = 1:K
            for j = 1:K
                for k = 1:K
                    temp2 = temp2 + prod(phi([i,j,k]))*T(i,j,k);
                end
            end
        end
        lambdaPerStart(tau) = temp2;
        phiPerStart(:,tau) = phi;
    end % end of random start tau
    % Find max component
    maxLambda = max(lambdaPerStart);
    taustar = find(lambdaPerStart==maxLambda, 1 );
    lambdaArr = [lambdaArr, maxLambda];
    phiMat = [phiMat, phiPerStart(:,taustar)];
    avgPowerIter = mean(powerIterCounts);
end

disp('Tensor power method complete!');

end