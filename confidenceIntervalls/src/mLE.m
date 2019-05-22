function [beta, fval] = mLE(dataX, dataY, modelfun, x0, additionalConstraints)
    params = size(x0,1);
    I1 = 1:params;
    I2 = (params+1):(params*2);
    I3 = (2*params+1):(params*3);
    lb = [zeros(params,1); -1000*ones(params,1); zeros(params,1)];
    ub = 1000*ones(3*params,1);

    %AddAdditional constraints;
    [Aeq, beq] = initConstraint(additionalConstraints);
    filter = not(isnan(additionalConstraints));
    x0(filter) = additionalConstraints(filter);
    
    betaEstimate = [x0(:,1); x0(:,2); x0(:,3)];
    
    paramFunction = @(beta) minFunction(beta, dataX, dataY, modelfun);
    
    %referenceValue = paramFunction(betaEstimate);
    
    [beta, fval] = getOptimalValue(paramFunction, betaEstimate, Aeq, beq, lb, ub);
    
    beta = [beta(I1) beta(I2) beta(I3)];
end
