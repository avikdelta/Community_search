function [thetaMatHat, runtime] = SpectralClusterNg(A,K,NumIter,fplot)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Implementation of spectral clustering algorithm of Ng, Jordan and Weiss.
% Implemented by Avik Ray (avik@utexas.edu)
% 
% Input parameters:
% A -- (n x n) adjacency matrix of graph G
% K -- number of clusters
% NumIter -- maximum number of iterations for k-mean clustering
% fplot -- flag for displaying plots
% 
% Output parameters:
% thetaMatHat -- (n x K) cluster membership matrix with 0/1 values
% runtime -- runtime in seconds
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic;

n = size(A,1);

% Compute degree
degree = sum(A,2);
D = diag(degree);

%% Compute Laplacian & svd

%L = D^(-.5)*(D-A)*D^(-.5);
%[U Sigma V] = svd(L);
%T = U(:,end-K+1:end);


L = D^(-.5)*(A+eye(n))*D^(-.5);
[U Sigma V] = svds(L,K);
T=U;


% Normalize rows
for i = 1:n
    T(i,:) = T(i,:)/norm(T(i,:));
end

%% k-mean cluster
%[M] = KMeans(T,K,NumIter,fplot);
%[M] = kmeans(T,K);

%% SLINK cluster
Z = linkage(T,'single','euclidean');
M = cluster(Z,'maxclust',K);

%% Assign clusters
thetaMatHat = zeros(n,K);
for k = 1:K
    kidx = find(M==k);
    thetaMatHat(kidx,k) = ones(length(kidx),1);
end

runtime = toc;

%% Plot community membership
if fplot
    figure; stem(M); grid on;
end

end
