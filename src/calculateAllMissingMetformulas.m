function model = calculateAllMissingMetformulas(model)
    [balanceStructure.elements, useMat, exitFlag]=parseFormulas(model.metFormulas, true);
    elementCode = balanceStructure.elements.abbrevs;
    
    missingFormulas = find(cellfun(@isempty,model.metFormulas));
    newMissing = length(missingFormulas);
    nrMissing = newMissing + 1;
    
    fprintf('Nr of missing formulas: %2.0f\n\n', newMissing);
    
    %iterate untill no further improvement
    while newMissing<nrMissing
        nrMissing = newMissing;
        affectedRxns = find(sum(abs(model.S(missingFormulas,:)),1));
        
        for i = 1:length(affectedRxns)
            model = tryToFixFormula(model, affectedRxns(i),missingFormulas, elementCode, useMat);
        end
        
        missingFormulas = find(cellfun(@isempty,model.metFormulas));
        [balanceStructure.elements, useMat, exitFlag]=parseFormulas(model.metFormulas, true);        
        newMissing = length(missingFormulas);
        fprintf('-----------\n');
    end
    fprintf('Nr of remaining missing formulas: %2.0f\n\n', newMissing);
end

function model = tryToFixFormula(model, rxn, missingFormulas, elementCode, useMat)
    includedMets = find(model.S(:,rxn));
    
    %We can fix the reaction if only one met is missing a formula
    target = intersect(includedMets,missingFormulas);
    
    if length(target) == 1
        stochiometry = full(model.S(includedMets,rxn));
        stochiometryMatrix = repmat(stochiometry, 1, length(elementCode));
        elementalBalance = useMat(includedMets, :).*stochiometryMatrix;
        delta = abs(sum(elementalBalance,1));
        formulaStr = makeFormulaString(delta, elementCode);
        if not(isempty(formulaStr))
            model.metFormulas{target} = formulaStr;
            fprintf('%s\t%s\n', model.metNames{target}, formulaStr)
        end
    end
end


function frmlStr = makeFormulaString(delta, elementCode)
    frmlStr = '';
    for i = 1:length(delta)
        if delta(i)>0
            frmlStr = [frmlStr elementCode{i} num2str(delta(i))];
        end
    end
end

    
    