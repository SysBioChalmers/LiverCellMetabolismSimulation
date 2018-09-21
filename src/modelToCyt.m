function modelToCyt(model, fluxes)
    model = removeCompartments(model);
    
    banRxns = {'PhosphatidylPool', 'humanGrowhOut', 'human_ATPMaintainance', 'human_DNAPool', 'human_RNAPool', 'human_proteinPool', 'lipidPool'};
    model.S(:,ismember(model.rxns, banRxns)) = 0;
    
    conLim = 6;
    
    conectivity = sum(abs(sign(model.S)),2);
    model.S(conectivity>conLim,:) = 0;
    

    
    adjMatrix = model.S * model.S';
    adjMatrix(eye(length(adjMatrix))==1) = 0;
    
    metNames = model.metNames;
    for i = 1:length(metNames)
        connectedMets = find(adjMatrix(i,:));
        for j=1:length(connectedMets)
            fprintf('%s\t%s\n', metNames{i}, metNames{connectedMets(j)});
        end
        adjMatrix(i,:) = 0;
        adjMatrix(:,i) = 0;
    end
end

function model = removeCompartments(model)
    metabolites = unique(model.metNames);
    for i = 1:length(metabolites)
        curMet = metabolites{i};
        matchingMets = find(ismember(model.metNames, curMet));
        model.S(matchingMets(1),:) = sum(model.S(matchingMets,:),1);
        model.S(matchingMets(2:end),:) = 0;
    end
end

