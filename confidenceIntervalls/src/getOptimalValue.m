function [beta, fval] = getOptimalValue(paramFunction, betaInValue, Aeq, beq, lb, ub)
    %'Display','off',
    options = optimoptions('fmincon', 'Display','off', 'MaxFunctionEvaluations', 1000000, 'MaxIterations', 1000000, 'OptimalityTolerance', 10^-10);
    [beta, fval] = fmincon(paramFunction, betaInValue, [], [], Aeq, beq, lb, ub, [], options);
end