function [Theta] = GenNonOverlapCluster(n,K,type,fsave,sizeArr,key)
% Script to generate community membership file
% [Theta] = GenNonOverlapCluster(1000,10,1,1);

Theta = zeros(n,K);
disp('Generating clusters ...');
if type == 1
    % All equal cluster
    typeId = 'Eq';
    size = floor(n/K);
    idx = 0;
    for i = 1:n
        if i > idx*size
            idx = idx + 1;
        end
        Theta(i,idx) = 1;
    end
elseif type == 2
    % Custom sizes
    typeId = 'Neq';
    startIdx = 1;
    for k = 1:K
        endIdx = startIdx + sizeArr(k) - 1;
        Theta(startIdx:endIdx,k) = ones(sizeArr(k),1);
        startIdx = endIdx + 1;
    end
end

if fsave
    fname = strcat('NOComm_', typeId, '_n', num2str(n),'_K', num2str(K), '_', key);
    comm.Theta = Theta;
    comm.n = n;
    comm.K = K;
    comm.sizes = sum(Theta,1);
    comm.type = type;
    comm.typeId = typeId;
    comm.sizeArr = sizeArr;
    comm.key = key;
    save(fname,'comm');
    disp('Community File saved!');
end

end