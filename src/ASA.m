function [normalizedGrowth, xValues, reactionNumbers] = ASA(model, fluxes, simSteps, reversible, rxnNrs, exceptions)
% Acute Sensitivity Analysis
% simulates the acute effect of pertubation of reactions
%
% [normalizedGrowth, xValues] = ASA(model, fluxes, simSteps)
%
% model                 a model in raven format
% fluxes                the flux vector with reference condition
% simSteps              the number of flux reduction steps to try 
% reversible            boolean, may the flux direction be reversed for reversible reactions
%---
% normalizedGrowth      sensitivity matrix, rows = rxns, cols = sensitivity
% xValues               flux reduction levels
% reactionNumbers       the numbers of the investigated fluxes
%
%
% Avlant Nilsson 2018
%
if ~exist('rxnNrs','var')
    rxnNrs=[];
end

if ~exist('exceptions','var')
    exceptions=[];
end

exceptionVals = fluxes(exceptions);

tresh = 10^-6;
epsilon = 10^-5;
growthTolerance = 10^-6;
preFilter = false;
dontRelaxMaintainance = true;

[exchangeRxns, exchangeRxnsIndexes]=getExchangeRxns(model,'both');

%%%%%%%%%%%%%%%%%
%Remove reactions that we are not interested in analyzing
if isempty(rxnNrs)
    toAnalyse = abs(fluxes)>tresh; %only reactions with flux
else
    toAnalyse = zeros(length(fluxes),1);
    toAnalyse(rxnNrs) = 1;
end


%Fix objective function 
objRxn = find(model.c);


if preFilter
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Find the obvious unconstrained candidates in parallell to save time
    %Minimize all affected fluxes
    model = setParam(model, 'lb', objRxn, fluxes(objRxn)*(1-growthTolerance));
    model = setParam(model, 'ub', objRxn, fluxes(objRxn));

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

    model.c = zeros(length(fluxes),1);
    model.c(objRxn) = 1;
end

%We should now be able to relax exchange bounds (not protein limited)
model.lb(exchangeRxnsIndexes) = -1000;
model.ub(exchangeRxnsIndexes) = 1000;

if dontRelaxMaintainance
    %Store maintainance information
    maintrxn = find(model.lb>0); %maintainance has a positive lower bond
    exceptions = [exceptions; maintrxn];
end

exceptionLB = model.lb(exceptions);
exceptionUB = model.ub(exceptions);

%%%%%%%%%%%%%%%%
%Prevent metabolic reorganization asside from utilized metabolic fluxes 
%i.e. prevents adaptation of proteom.
possitiveRxns = fluxes>tresh;
reversibleReactions = model.lb<0;

model.ub(possitiveRxns) = fluxes(possitiveRxns);
model.lb(possitiveRxns) = 0;
if reversible
    model.lb(and(possitiveRxns, reversibleReactions)) = -fluxes(and(possitiveRxns, reversibleReactions));
end

neutralRxns = abs(fluxes)<tresh;
model.ub(neutralRxns) = 0;
model.lb(neutralRxns) = 0;    

negativeRxns = fluxes<-tresh;
forwardReactions = model.ub>0;
model.lb(negativeRxns) = fluxes(negativeRxns);
model.ub(negativeRxns) = 0; 
if reversible
    model.ub(and(negativeRxns, forwardReactions)) = -fluxes(and(negativeRxns, forwardReactions));   
end
%Reset exceptions if any
model.lb(exceptions) = exceptionLB;
model.ub(exceptions) = exceptionUB;
    
%We should now be able to relax exchange bounds (not protein limited)
model.lb(exchangeRxnsIndexes) = -1000;
model.ub(exchangeRxnsIndexes) = 1000;

%Sensitivity analysis
reactionNumbers = find(toAnalyse);
model = setParam(model, 'lb', objRxn, 0);
model = setParam(model, 'ub', objRxn, 1000);

simulationSteps = linspace(0,1, simSteps);
simulationSteps(end) = []; %we already know reference value

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
        
        if curAmount < 0
            model = setParam(model, 'lb', curRxn, curAmount*(1+epsilon));
            model = setParam(model, 'ub', curRxn, curAmount*(1-epsilon));
        else
            model = setParam(model, 'lb', curRxn, curAmount*(1-epsilon));
            model = setParam(model, 'ub', curRxn, curAmount*(1+epsilon));            
        end

        solution = solveLin(model,1);
        growthRates(i,j) = -solution.f;        

%         if isempty(solution.x)
%             curAmount
%             constructEquations(model, curRxn)
%            j 
%         end

%           if j == 2
%              solution = solveLinMin(model,1);
%              printFluxesAvlant(model, 1000*solution.x, true) 
%           end
    end
    %Reset bounds
   model.ub = allUbs;
   model.lb = allLbs;
end



normalizedGrowth = growthRates/referenceGrowth;

%no restriction obviously gives maximum growth
normalizedGrowth = [normalizedGrowth ones(length(reactionNumbers),1)];
xValues = [simulationSteps 1];

end

