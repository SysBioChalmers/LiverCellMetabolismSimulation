function [conf] = mLEConfidence(dataX, dataY, modelfun, beta, priorSTD)
    params = size(beta,1);
    betaIn = [beta(:,1); beta(:,2); beta(:,3)];
    fastCI = true;

    confidencebound = 0.341 * 2;
    
    conf = zeros(params,2);
    
    %unconstrained value
    fvalRef = minFunction(betaIn, dataX, dataY, modelfun);

    paramFunction = @(beta) minFunction(beta, dataX, dataY, modelfun);
    
    %Add Prior STD constraints
    Aeq = zeros(1, length(betaIn));
    beq = 0;
    for i = 1:length(priorSTD)
        curPar = params*2 + i;
        if priorSTD(i)~=0
            [Aeq, beq] = addEqConstraints(Aeq, beq, length(betaIn), curPar, priorSTD(i));
            betaIn(curPar) = priorSTD(i);
        end
    end
    
    %Calculate CI for growth using all available data
    conf(1,:) = getCIofVariable(paramFunction, betaIn, Aeq, beq, params + 1, fvalRef, confidencebound);
    

    for i = 2:params
        disp(i)
        if fastCI
            %Calculate CI for metabolites using only data for growth and the specific amino acid
            %This generates around larger CI due to missing interaction effects
            %with growth, but is many times faster.      
            conf(i,:) = calculateFastCI(dataX, dataY, modelfun, betaIn, Aeq, beq, params + i, confidencebound);
        else
            conf(i,:) = getCIofVariable(paramFunction, betaIn, Aeq, beq, params + i, fvalRef, confidencebound);
        end    
    end    
end

function CI = calculateFastCI(dataX, dataY, modelfun, betaIn, Aeq, beq, parameter, confidencebound)
    params = size(betaIn,1)/3;
    curData = parameter - params;
    valueOfAAData = 2;
    paramLength = 2;
    
    %Remove independent data
    filter = ismember(dataX(:,2), [1, curData]);
    curDataX = dataX(filter,:);
    curDataY = dataY(filter,:);
    curDataX(curDataX(:,2)==curData,2) = valueOfAAData;    

    %Remove independet variables
    varablesToKeep = [(1) (curData) (params+1) (params+curData) (params*2+1) (params*2+curData)];
    betaIn = betaIn(varablesToKeep,:);
    Aeq = Aeq(:,varablesToKeep);
    %Beq would also need adjusting if Aeq contained any values not in
    %growth or curData
    
    fvalRef = minFunction(betaIn, curDataX, curDataY, modelfun);
    paramFunction = @(beta) minFunction(beta, curDataX, curDataY, modelfun);

    CI = getCIofVariable(paramFunction, betaIn, Aeq, beq, paramLength+valueOfAAData, fvalRef, confidencebound);
    if length(CI)<2
       CI = zeros(:,2); 
    end
end


function CI = getCIofVariable(paramFunction, betaIn, Aeq, beq, parameter, fvalRef, confidencebound)
    global betaInValue %to update staringpoint of solver between each iteration
    betaInValue = betaIn;
    refParameter = betaIn(parameter);
    
    fvalTresh = fvalRef + chi2inv(confidencebound, 1);
    
    %Fix parameter for effect investigation and ensure feasibility of startingpoint
    [Aeq, beq] = addEqConstraints(Aeq, beq, length(betaIn), parameter, refParameter);
    
    
    LB = true;
    searchProblem = @(x) testConstrainedParam(x, LB, fvalTresh, paramFunction, parameter, refParameter, Aeq, beq);
    x0 = 0.8 * refParameter;
    [lb, f] = fminsearch(searchProblem, x0);
    
    LB = false;
    searchProblem = @(x) testConstrainedParam(x, LB, fvalTresh, paramFunction, parameter, refParameter, Aeq, beq);
    x0 = 1.2 * refParameter;
    [ub, f] = fminsearch(searchProblem,x0);

    CI = [lb ub];
end

function fitness = testConstrainedParam(curVal, LB, fvalTresh, paramFunction, parameter, refParameter, Aeq, beq)
    global betaInValue
    options = optimoptions('fmincon', 'Display','off', 'MaxFunctionEvaluations', 100000, 'MaxIterations', 100000, 'OptimalityTolerance', 10^-8);
    
    penaltyForWrongDirection = 0;
    if and(LB == true, refParameter<curVal)
        penaltyForWrongDirection = (refParameter-curVal)^2;
    elseif and(LB == false, refParameter>curVal)
        penaltyForWrongDirection = (refParameter-curVal)^2;
    end
    
    beq(length(beq)) = curVal; %assumes latest changed beq
    betaInValue(parameter) = curVal;
    
    [betaInValue,fval] = fmincon(paramFunction, betaInValue, [], [], Aeq, beq, [], [], [], options);
    
    fitness = (fval-fvalTresh)^2 + penaltyForWrongDirection; 
end

function [Aeq, beq] = addEqConstraints(Aeq, beq, nrOfVariables, variable, value)
    newAeq = zeros(1,nrOfVariables);
    newAeq(variable) = 1;
    Aeq = [Aeq;newAeq];
    beq = [beq;value];
end

function plotCI(xvalues, mLEResponse, fvalRef, fvalTresh, xi)
    %xvalues = linspace(0.3,1.8, 30);
    %mLEResponse = zeros(length(xvalues),1);
    %for i = 1:length(xvalues)        
    %plotCI(xvalues, mLEResponse, fvalRef, fvalTresh, xi);    
    %end
    %xi = polyxpoly(xvalues,mLEResponse,xvalues,ones(1,length(xvalues))*fvalTresh);

    hold all
    plot(xvalues, mLEResponse)
    plot(1, fvalRef, 'ko')
    plot([0.3 1.8],  fvalRef * [1 1], 'k--')
    plot([0.3 1.8],  fvalTresh * [1 1], 'k--')
    plot(xi, fvalTresh * ones(length(xi),1), 'ko')  
end