function [A] = GenSBMGraph(Theta,p,q,fplot)

[n K] = size(Theta);

% Generate edges
disp('Generating graph ...');
A = zeros(n);
for i = 1:n
    for j = 1:i-1
        if Theta(i,:)*Theta(j,:)' == 0
            A(i,j) = double(rand<q);
        else
            A(i,j) = double(rand<p);
        end
        A(j,i) = A(i,j);
    end
end

% Display stats and plots
if fplot
    % Compute degree
    Deg = sum(A,1);
    maxdeg = max(Deg);
    mindeg = min(Deg);

    disp(['Max degree = ' num2str(maxdeg) ' Min degee = ' num2str(mindeg)]);
    NumberMinDeg = length(find(Deg==mindeg));
    disp(['Number of nodes with minimum degree = ' num2str(NumberMinDeg)]);
    hist(Deg,1:maxdeg);
    
    % Compute spectrum
    L = diag(Deg) - A;
    A_dash = A + eye(n);
    [U,S,V] = svd(A_dash);
    figure; plot([1:n],diag(S),'r-o'); grid on;
    minSVal = min(diag(S));
    disp(['Minimum sigular value = ' num2str(minSVal)]);
    
    c_node_count = sum(Theta,1);
    nEdges = sum(Deg)/2;
    disp(['Number of edges = ' num2str(nEdges)]);
end

end