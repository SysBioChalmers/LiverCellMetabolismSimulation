function model = bindFBA(model, mets, fluxes)
    [id, exchangeRxns] =getExchangeRxns(model);
    reactionNumbers = getBounds(model, mets);
    model = setParam(model, 'lb', exchangeRxns, 0);
    model = setParam(model, 'ub', exchangeRxns, 1000);    
    model = setParam(model, 'lb', reactionNumbers, fluxes);
    model = setParam(model, 'ub', reactionNumbers, fluxes);
    
    provideCholine = true;
    provideUnused = false;
    provideInositol = true;
    
    upMets = {'pantothenate[s]', 'nicotinamide[s]', 'H2O[s]', 'O2[s]', 'H+[s]', 'NH3[s]', 'urea[s]', 'Pi[s]', 'Fe2+[s]', 'sulfate[s]', 'taurine[s]'};
    uptake = getBounds(model, upMets);
    model = setParam(model, 'lb', uptake, -1000);

    reMets = {'HCO3-[s]', 'CO2[s]', 'H+[s]', 'NH3[s]', 'urea[s]', 'H2S[s]', 'biomass[s]'};
    release = getBounds(model, reMets);
    model = setParam(model, 'ub', release, 1000);
    
    objectiveFunction = 'HumanGrowth';
    model = setParam(model, 'obj', objectiveFunction, 1);
    model = setParam(model, 'lb', objectiveFunction, 0);
    model = setParam(model, 'ub', objectiveFunction, 1000);


    
    if provideCholine
        %block choline catabolism
        model = setParam(model, 'ub', 'HMR_4696', 0);
        model = setParam(model, 'ub', 'HMR_8441', 0);
        upMets = {'choline[s]'};
        uptake = getBounds(model, upMets);
        model = setParam(model, 'lb', uptake, -1000);
    end
    
    if provideInositol
       %block inositol catabolism
        model = setParam(model, 'ub', 'HMR_6539', 0);
        upMets = {'inositol[s]'};
        uptake = getBounds(model, upMets);
        model = setParam(model, 'lb', uptake, -1000);        
    end
    
    if provideUnused
        upMets = {'riboflavin[s]', 'folate[s]'};
        uptake = getBounds(model, upMets);
        model = setParam(model, 'lb', uptake, -1000);           
    end
    
end

