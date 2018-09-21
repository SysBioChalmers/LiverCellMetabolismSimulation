function solution = minimizeFluxThroughReactions(model, reactions)
    nrOfReactions = length(model.c);   
       
    %Add reversed Matrix
    newModel = model;
    newModel.S = horzcat(model.S, -model.S);    
    newModel.lb = vertcat(model.lb,  -model.ub);
    newModel.ub = vertcat(model.ub, -model.lb);
    newModel.c = vertcat(model.c, -model.c);

    %Set Lower Bound
    newModel.lb(newModel.lb<0) = 0;
    newModel.ub(newModel.ub<0) = 0;

    %Minimize fluxes
    newModel.c = zeros(nrOfReactions*2,1);
    newModel.c(reactions) = -1;
    newModel.c(reactions + nrOfReactions) = -1;
    
    solution = solveLin(newModel, true);
    solution = solution.x(1:nrOfReactions) - solution.x((nrOfReactions+1):end);
end