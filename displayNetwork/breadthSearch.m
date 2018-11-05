function [l, sp] = breadthSearch(A,s,d)
maxDepth = 6;

l = inf;
sp = [];

if s==d
   l=1;
   sp = [s];
else
    for i = 2:maxDepth
        sVector = zeros(1, i);
        sVector(1) = s;
        paths = recursiveSearch(A, sVector, d);
        if not(isempty(paths))
            l = i;
            sp = paths;
            break
        end
    end
end

end

function sp = recursiveSearch(A, sVector, d)
    currentDepth = sum(sVector>0);
    sp = [];
    
    if currentDepth>=length(sVector)
        if sVector(end) == d
           sp = sVector; %Found d! 
        end
    else
        currentNode = sVector(currentDepth);
        conectingPaths = find(A(currentNode,:));
        tmpVector = sVector;
        for i = 1:length(conectingPaths)
            tmpVector(currentDepth+1) = conectingPaths(i);
            curPath = recursiveSearch(A, tmpVector, d);
            sp = [sp; curPath];
        end
    end
end