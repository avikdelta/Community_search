function communitySearchAll(CommMembershipFile, p, q, numLabel)
%--------------------------------------------------------------
% Main script to perform community search in a SBM given side information
% as few labeled nodes per community
%
% Author: Avik Ray (avik@utexas.edu)
%
% Input params:
% CommMembershipFile: mat file defining community membership structure info 
% p: intra-community edge probability
% q: inter-community edge probability
% numLabel: Number of labeled nodes (side-information) / community
%
% Example:
% communitySearchAll('NOComm_Eq_n1000_K10.mat', 0.1, 0.01, 10)
%--------------------------------------------------------------

close all hidden;

% Load community membership matrix
load(CommMembershipFile);

% initialize variables
Theta = comm.Theta;
[n K] = size(Theta);
commSizes = comm.sizes;
commStartArr = 1 + [0 cumsum(commSizes(1:K-1))];

% Generate stochastic block model graph
[A] = GenSBMGraph(Theta,p,q,0);
    
% Iterate over all community to search for each of them
muMatS = zeros(n,K);
if numLabel > 0
    labelNodeArr = [];
    clusterArr = [];
end
        
for k = 1:K
    
    % wlog, we can assume the first "numLabel" nodes in the community are labeled
    startidx = commStartArr(k);
    endidx = startidx + numLabel - 1;
    
    % initialize side-information vector as the average of one hop
    % neighbors adjacency vector
    v = sum(A(:,startidx:endidx),2)/numLabel;
    
    % store the labeled nodes
    labelNodeArr = [labelNodeArr startidx:endidx];
    clusterArr = [clusterArr k*ones(1,numLabel)];
            
    % Run community search
    disp(['Running community search algo for k = ' num2str(k)]);
    tau = 0; % set decision threshold as 0, since we will perform a greedy clustering
    [thetaHat, inrprod, runtime] = searchSBMWhiten(A,v,K,p,q,tau,0);
     muMatS(:,k) = inrprod;
end

% Perform greedy clustering of nodes
thetaMatW = zeros(n,K);
for i = 1:n
    lidx = find(labelNodeArr==i);
    if ~isempty(lidx)
        % labeled node
        kidx = clusterArr(lidx);
        thetaMatW(i,kidx) = 1;
    else
        % assign node to the community with which the adjacency vector has
        % the highest inner product
        proj = A(i,:)*muMatS;
        kidx = find(proj==max(proj),1);  
        thetaMatW(i,kidx) = 1; 
     end
end

% Evaluate error
errVec = zeros(K, 1);
for k = 1:K
    % note that the columns of thetaMatW can be permuted
    [error permidx] = ClusterErr(Theta, thetaMatW, k);
    errVec(k, 1) = error;
end

total_error = sum(errVec)/n;
disp(['Percentage error = ', num2str(total_error*100), '%']);

disp('Experiment complete !');

end