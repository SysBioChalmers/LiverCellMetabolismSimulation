function [metNames, percents] = printAnaplerosis(smallModel, smallSolution)
    metRxn = find(contains(smallModel.rxns,'metabolitePool'));
    poolMets = find(smallModel.S(:,metRxn));
    poolFlux = smallSolution(metRxn);

    metNames = smallModel.metNames(poolMets);
    percents = zeros(length(poolMets),1);
    for i = 1:length(poolMets)
       curMet = poolMets(i);
       metName = smallModel.metNames{curMet};
       curRxns = find(smallModel.S(curMet,:)>0);
       curFlux = abs(poolFlux * full(smallModel.S(curMet, metRxn)));
       curStochiometry = full(smallModel.S(curMet,curRxns));
       totFlux = sum(curStochiometry' .* smallSolution(curRxns));
       percents(i) = curFlux./totFlux;
       fprintf('%s\t%2.6f\t%2.6f\t%2.6f\n', metNames{i}, curFlux, totFlux, percents(i));
    end
end

