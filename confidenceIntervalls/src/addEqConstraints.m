function [Aeq, beq] = addEqConstraints(Aeq, beq, nrOfVariables, variable, value)
    newAeq = zeros(1,nrOfVariables);
    newAeq(variable) = 1;
    Aeq = [Aeq;newAeq];
    beq = [beq;value];
end