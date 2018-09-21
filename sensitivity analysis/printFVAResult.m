function printFVAResult(model, result)
fprintf('rxn\teq\tsub\tlb\tflux\tub\tDiffSize\n');

    for i = 1:length(result)
        rxnId = result(i,1);
        rxn = model.rxns{rxnId};
        eq = constructEquations(model, rxnId);
        sub = model.subSystems{rxnId};    
        flux = 1000 * result(i,2);        
        lb = 1000 * result(i,3);
        ub = 1000 * result(i,4);
        percentDiff = 100*abs(ub-lb)/2/abs(flux);
        fprintf('%s\t%s\t%s\t%2.2f\t%2.2f\t%2.2f\t%2.2f\n', rxn, eq{1}, sub, lb, flux, ub, percentDiff);
    end
end

