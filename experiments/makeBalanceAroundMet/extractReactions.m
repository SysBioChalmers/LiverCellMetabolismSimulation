function [inRxns, outRxns, influx, outflux] = extractReactions(model, smallSolution, investigateMetabolite)
    treshold = 10^-6;
    metId = findIndex(model.metNames, investigateMetabolite);
    reactions = find(sum(abs(model.S(metId,:)),1));
    stochiometry = full(model.S(metId,reactions))';
    totalFlux = stochiometry .* smallSolution(reactions);
    totalFlux = sum(totalFlux, 2);

    inRxns = totalFlux>treshold;
    outRxns = totalFlux<-treshold;
  
    influx = totalFlux(inRxns);
    outflux = totalFlux(outRxns);
    
    influx = abs(influx);
    [influx, indx] = sort(influx, 'descend');
    inRxns = reactions(inRxns);
    inRxns = inRxns(indx);    
    
    outflux = abs(outflux);
    [outflux, indx] = sort(outflux, 'descend');
    outRxns = reactions(outRxns);
    outRxns = outRxns(indx);
    
    
end

