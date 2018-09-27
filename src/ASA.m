function [normalizedGrowth, xValues, reactionNumbers] = ASA(model, fluxes, simSteps)
% Acute Sensitivity Analysis
% simulates the acute effect of pertubation of reactions
%
% [normalizedGrowth, xValues] = ASA(model, fluxes, simSteps)
%
% model                 a model in raven format
% fluxes                the flux vector with reference condition
% simSteps              the number of flux reduction steps to try 
% --
% normalizedGrowth      sensitivity matrix, rows = rxns, cols = sensitivity
% xValues               flux reduction levels
% reactionNumbers       the numbers of the investigated fluxes
%
% Avlant Nilsson 2018
%

tresh = 10^-6;
growthTolerance = 0;
preFilter = false;
dontRelaxMaintainance = true;

[exchangeRxns, exchangeRxnsIndexes]=getExchangeRxns(model,'both');

%%%%%%%%%%%%%%%%%
%Remove reactions that we are not interested in analyzing
toAnalyse = abs(fluxes)>tresh; %only reactions with flux

%Fix objective function 
objRxn = find(model.c);
model = setParam(model, 'lb', objRxn, fluxes(objRxn)*(1-growthTolerance));
model = setParam(model, 'ub', objRxn, fluxes(objRxn));

if preFilter
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Find the obvious unconstrained candidates in parallell to save time
    %Minimize all affected fluxes
    model.c = zeros(length(fluxes),1);
    model.c(toAnalyse) = -1;
    curSol = solveLin(model,1);
    lowerBound = abs(curSol.x)>900;

    %Maximize all fluxes
    model.c = zeros(length(fluxes),1);
    model.c(toAnalyse) = 1;
    curSol = solveLin(model,1);
    upperBound = abs(curSol.x)>900;

    toAnalyse(and(lowerBound, upperBound)) = false;
end



if dontRelaxMaintainance
    %Store maintainance information
    exceptions = model.lb>0; %maintainance has a positive lower bond
    exceptionVals = model.lb(exceptions);
else
    exceptions = [];
    exceptionVals = [];
end

%%%%%%%%%%%%%%%%
%Prevent metabolic reorganization asside from utilized metabolic fluxes 
%i.e. prevents adaptation of proteom.
possitiveRxns = sign(fluxes)>0;
model.ub(possitiveRxns) = fluxes(possitiveRxns);
model.lb(possitiveRxns) = 0;

neutralRxns = sign(fluxes)==0;
model.ub(neutralRxns) = 0;
model.lb(neutralRxns) = 0;    

negativeRxns = sign(fluxes)<0;
model.lb(negativeRxns) = fluxes(negativeRxns);
model.ub(negativeRxns) = 0;    

%Reset exceptions if any
model.lb(exceptions) = exceptionVals;

%We should be able to relax exchange bonds as well (not protein limited)
model.lb(exchangeRxnsIndexes) = -1000;
model.ub(exchangeRxnsIndexes) = 1000;
    

%Sensitivity analysis
reactionNumbers = find(toAnalyse);
model.c = zeros(length(fluxes),1);
model.c(objRxn) = 1;
model = setParam(model, 'lb', objRxn, 0);
model = setParam(model, 'ub', objRxn, 1000);

simulationSteps = linspace(0,1, simSteps);
simulationSteps(end) = []; %we already know default value

%Reference
referenceAmounts = fluxes(reactionNumbers);
referenceGrowth = model.c' * fluxes;
allUbs = model.ub;
allLbs = model.lb;

growthRates = zeros(length(reactionNumbers),length(simulationSteps));

for i = 1:length(reactionNumbers)
    i
    curRxn = reactionNumbers(i);    
    for j = 1:length(simulationSteps)
        curAmount = referenceAmounts(i) * simulationSteps(j);
        model = setParam(model, 'lb', curRxn, curAmount);
        model = setParam(model, 'ub', curRxn, curAmount);
        solution = solveLin(model,1);
        growthRates(i,j) = -solution.f;
    end
    
    %Reset bounds
   model.ub = allUbs;
   model.lb = allLbs;
end

normalizedGrowth = growthRates/referenceGrowth;

%Obviously no restriction gives maximum growth
normalizedGrowth = [normalizedGrowth ones(length(reactionNumbers),1)];
xValues = [simulationSteps 1];

end

