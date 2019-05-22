function [Aeq, beq] = initConstraint(aditionalConstraints)
    addConstr = [aditionalConstraints(:,1); aditionalConstraints(:,2); aditionalConstraints(:,3)];
    Aeq = [];
    beq = []; 
    
    for i = 1:length(addConstr)
        if not(isnan(addConstr(i)))
            [Aeq,beq] = addEqConstraints(Aeq, beq, length(addConstr), i, addConstr(i));
        end
    end
end

