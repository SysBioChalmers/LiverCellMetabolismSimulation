function candidates = getGrowthLimitingMet(model, fluxMets, fluxValues)
    epsilon = 0.0001;
    tolerance = 10^-6;

    growthRate = getGrowthRate(model, fluxMets, fluxValues);
    
    
    candidates = 1:length(fluxMets);
    
    while not(isempty(candidates))
        midPoint = ceil(length(candidates)/2);
        newCandidates = candidates(1:midPoint);
        testFluxes = fluxValues;
        testFluxes(newCandidates) = testFluxes(newCandidates) * (1-epsilon);
        testModel = model;
        testModel = bindFBA(testModel, fluxMets, testFluxes);
        testGrowth = getGrowthRate(testModel, fluxMets, testFluxes);
        
        
        if (testGrowth/growthRate) < (1-epsilon)+tolerance
            candidates = newCandidates;
            if length(candidates) == 1
                break
            end
        else
            if length(candidates) == 1
                candidates = [];
                break
            else
                candidates = candidates((midPoint+1):end);
            end       
        end

    end
    
fluxMets(candidates)    

end

function growthRate = getGrowthRate(model, fluxMets, fluxValues)
    model = bindFBA(model, fluxMets, fluxValues);
    growthRate = runFBA(model, {'biomass[s]'}, true);
end