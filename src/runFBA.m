function fluxes = runFBA(model, mets, verbose)
    objectiveFunction = 'HumanGrowth';
    model = setParam(model, 'obj', objectiveFunction, 1);
    solution = solveLinMin(model, not(verbose));
    
    reactionNumbers = getBounds(model, mets);
    if length(solution.x) == 1
        if verbose
            disp('warning no solution found')
        end
        fluxes = zeros(length(reactionNumbers),1);
    else
        fluxes = solution.x(reactionNumbers);
    end
end

