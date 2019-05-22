function [betaRef, conf] = mLEConfidence(dataX, dataY, modelfun, beta, confidencebound, additionalConstraints)
    params = size(beta,1);
    I1 = 1:params;
    I2 = (params+1):(params*2);
    I3 = (2*params+1):(params*3);
    lb = [zeros(params,1); -1000*ones(params,1); zeros(params,1)];
    ub = 1000*ones(3*params,1);    
    
    betaIn = [beta(:,1); beta(:,2); beta(:,3)];
    
    
    conf = zeros(params,2);
    
    %AddAdditional constraints;
    [Aeq, beq] = initConstraint(additionalConstraints);
    
    %get reference value
    paramFunction = @(beta) minFunction(beta, dataX, dataY, modelfun);
    [betaRef, fvalRef] = getOptimalValue(paramFunction, betaIn, Aeq, beq, lb, ub);
       
    %Calculate CI for growth
    conf(1,:) = getCIofVariable(paramFunction, betaRef,  params + 1, fvalRef, confidencebound, Aeq, beq, lb, ub);
    
    for i = 2:params
        conf(i,:) = getCIofVariable(paramFunction, betaRef, params + i, fvalRef, confidencebound, Aeq, beq, lb, ub); 
    end
    
    betaRef = [betaRef(I1) betaRef(I2) betaRef(I3)];
end

function CI = getCIofVariable(paramFunction, betaIn, parameter, fvalRef, confidencebound, Aeq, beq, lb, ub)
    global betaInValue %to store staringpoint of solver between each iteration
    betaInValue = betaIn;
    refParameter = betaIn(parameter);
    options = optimset('MaxFunEvals',1000, 'MaxIter', 2000);
 
    fvalTresh = fvalRef + chi2inv(confidencebound, 1);
    
    %Fix parameter for effect investigation and ensure feasibility of startingpoint
    [Aeq, beq] = addEqConstraints(Aeq, beq, length(betaIn), parameter, refParameter);
    
    LBsetting = true;
    searchProblem = @(x) testConstrainedParam(x, LBsetting, fvalTresh, paramFunction, parameter, refParameter, Aeq, beq, lb, ub);
    x0 = 0.8 * refParameter;
    [CIlb, f, exitflag] = fminsearch(searchProblem, x0, options);
    
    LBsetting = false;
    searchProblem = @(x) testConstrainedParam(x, LBsetting, fvalTresh, paramFunction, parameter, refParameter, Aeq, beq, lb, ub);
    x0 = 1.2 * refParameter;
    [CIub, f, exitflag] = fminsearch(searchProblem, x0, options);
    
    CI = [CIlb CIub];
end

function fitness = testConstrainedParam(curVal, LB, fvalTresh, paramFunction, parameter, refParameter, Aeq, beq, lb, ub)
    global betaInValue
   
    penaltyForWrongDirection = 0;
    if and(LB == true, refParameter<curVal)
        penaltyForWrongDirection = (refParameter-curVal)^2;
    elseif and(LB == false, refParameter>curVal)
        penaltyForWrongDirection = (refParameter-curVal)^2;
    end
    
    beq(length(beq)) = curVal; %assumes latest changed beq
    betaInValue(parameter) = curVal;
    
    [betaInValue, fval] = getOptimalValue(paramFunction, betaInValue, Aeq, beq, lb, ub);
    
    fitness = (fval-fvalTresh)^2 + penaltyForWrongDirection; 
end




