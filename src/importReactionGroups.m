function [groupNames, reactionGroups, rxnStochiometry, coordinates] = importReactionGroups(fileName)
    raw = IO(fileName);
    groupNames = raw(2:end,1);
    for i = 1:length(groupNames)
        reactionGroups{i} = strsplit(raw{i+1,2},';');
    end

    if size(raw,2)>2
        for i = 1:length(groupNames)
            rxnStochiometry{i} = strsplit(raw{i+1,3},';');
            rxnStochiometry{i} = cell2nummat(rxnStochiometry{i});
        end    
    end

    if size(raw,2)>3
        for i = 1:length(groupNames)
            coordinates{i} = strsplit(raw{i+1,4},';');
            coordinates{i} = cell2nummat(coordinates{i});
        end    
    end    
end

