function [beta, fval] = mLE(dataX, dataY, modelfun, x0, std0, priorSTD)
    params = length(std0);
    I1 = 1:params;
    I2 = (params+1):(params*2);
    I3 = (2*params+1):(params*3);

    %'Display', 'iter', , 'Algorithm', 'quasi-newton', 
    options = optimoptions('fmincon', 'MaxFunctionEvaluations', 100000, 'MaxIterations', 100000, 'OptimalityTolerance', 10^-8);
    
    betaEstimate = [x0(:,1); x0(:,2); std0];
    
    paramFunction = @(beta) minFunction(beta, dataX, dataY, modelfun);

    %Add Prior STD constraints
    Aeq = [];
    beq = [];
    for i = 1:length(priorSTD)
        if priorSTD(i)~=0
            curPar = params*2 + i;
            [Aeq, beq] = addEqConstraints(Aeq, beq, length(betaEstimate), curPar, priorSTD(i));
            %betaEstimate(curPar) = priorSTD(i);
        end
    end
    
    Aeq

    [beta,fval] = fminunc(paramFunction, betaEstimate, options);
    %[beta,fval] = fmincon(paramFunction, betaEstimate, [], [], Aeq, beq, [], [], [], options);
    
    beta = [beta(I1) beta(I2) beta(I3)];
end



function [Aeq, beq] = addEqConstraints(Aeq, beq, nrOfVariables, variable, value)
    newAeq = zeros(1,nrOfVariables);
    newAeq(variable) = 1;
    Aeq = [Aeq;newAeq];
    beq = [beq;value];
end

