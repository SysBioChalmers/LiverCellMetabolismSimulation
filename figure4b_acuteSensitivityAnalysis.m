clc
clf
load('model/genericHuman2')
addpath('src')
addpath('sensitivity analysis')

[fluxMets, fluxValues] = loadFluxes('fluxvalues', 'hepg2-6mm-.txt');
model = setupBiomass(model, 150, 0.5);
model = bindFBA(model, fluxMets, fluxValues(:,2)/1000);

curSol = solveLinMin(model,1);
fluxes = curSol.x;
growthTolerance = 0;
tresh = 10^-6;

%%%%%%%%%%%%%%%%%%%%
%Relax upper bound for uptake reactions to allow sub maximal flux which may
%be required to prevent linear pathways from beeing infeasible
[rxn, id] = getExchangeRxns(model);
for i = 1:length(id)
    if model.lb(id(i)) < 0
        model.ub(id(i)) = 1000;
    end
end


%%%%%%%%%%%%%%%%%
%Remove reactions that we are not interested in analyzing
toAnalyse = abs(fluxes)>tresh; %only reactions with flux

%Fix objective function 
objRxn = find(model.c);
model = setParam(model, 'lb', objRxn, fluxes(objRxn)*(1-growthTolerance));
model = setParam(model, 'ub', objRxn, fluxes(objRxn));

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

%%%%%%%%%%%%%%%%
%Prevent metabolic reorganization asside from utilized metabolic fluxes 
%i.e. prevents adaptation of proteom.
subMaxModel = true;
if subMaxModel
    exceptions = model.lb>0; %maintainance has a positive lower bond
    exceptionVals = model.lb(exceptions);
    
    possitiveRxns = sign(fluxes)>0;
    model.ub(possitiveRxns) = fluxes(possitiveRxns);
    model.lb(possitiveRxns) = 0;
    
    %zero reacions
    neutralRxns = sign(fluxes)==0;
    model.ub(neutralRxns) = 0;
    model.lb(neutralRxns) = 0;    
    
    negativeRxns = sign(fluxes)<0;
    model.lb(negativeRxns) = fluxes(negativeRxns);
    model.ub(negativeRxns) = 0;    
    
    model.lb(exceptions) = exceptionVals;
    
    %We should be able to relax exchange bonds as well
    model.lb(id) = -1000;
    model.ub(id) = 1000;
end

%Sensitivity analysis
reactionNumbers = find(toAnalyse);
model.c = zeros(length(fluxes),1);
model.c(objRxn) = 1;
model = setParam(model, 'lb', objRxn, 0);
model = setParam(model, 'ub', objRxn, 1000);

simulationSteps = linspace(0,1, 10);
simulationSteps(end) = [];

growthRates = zeros(length(reactionNumbers),length(simulationSteps));

%Reference
solution = solveLinMin(model,1);
referenceAmounts = solution.x(reactionNumbers);
referenceGrowth = -solution.f

glucoseLactate = getBounds(model, {'glucose[s]', 'L-lactate[s]'});
glucoseLactateratio = solution.x(glucoseLactate(2))/glucoseLactate(1);

allUbs = model.ub;
allLbs = model.lb;

for i = 1:length(reactionNumbers)
    i
    curRxn = reactionNumbers(i);    
    for j = 1:length(simulationSteps)
        curAmount = referenceAmounts(i) * simulationSteps(j);
        model = setParam(model, 'lb', curRxn, curAmount);
        model = setParam(model, 'ub', curRxn, curAmount);
%         if curRxn == glucoseLactate(1) %if glucose then reduce lactate as well
%             model = setParam(model, 'lb', glucoseLactate(2), curAmount*glucoseLactateratio);
%             model = setParam(model, 'ub', glucoseLactate(2), curAmount*glucoseLactateratio);      
%         end
        solution = solveLin(model,1);
        growthRates(i,j) = -solution.f;
    end
    
    %Reset bounds
   model.ub = allUbs;
   model.lb = allLbs;
end

save('sensitivity analysis/sensitivityProfilesMaxFlux.mat','growthRates','reactionNumbers', 'simulationSteps')

normalizedGrowth = growthRates/referenceGrowth;
%Obviously no restriction gives maximum growth
normalizedGrowth = [normalizedGrowth ones(length(reactionNumbers),1)];
xValues = [simulationSteps 1];

hold all
%Displace the lines verticly at random by 1% to prevent overlapping lines
randValues = 0.02 * rand(length(reactionNumbers), 1)-0.01; 
normalizedGrowthAndDisplacement = normalizedGrowth + repmat(randValues, 1, length(xValues));

for i = 1:size(normalizedGrowthAndDisplacement,1)
    plot(xValues, normalizedGrowthAndDisplacement(i,:), 'color', [0.3 0.3 0.3 0.1]);
end

colors = [
         0    0.4470    0.7410
    0.8500    0.3250    0.0980
    0.9290    0.6940    0.1250
    0.4940    0.1840    0.5560
    0.4660    0.6740    0.1880
    0.3010    0.7450    0.9330
    0.6350    0.0780    0.1840];

conditions = {'glucose[s]', 'glutamine[s]', 'O2[s]'};
interestingReactions = getBounds(model, conditions);

%plot([0 1], [0, 1], '--', 'color', [0 0 0], 'linewidth', 2);


legendIds = zeros(length(interestingReactions),1);
for i = 1:length(interestingReactions)
    curRxn = ismember(reactionNumbers, interestingReactions(i));
    legendIds(i) = plot(xValues, normalizedGrowth(curRxn,:), 'linewidth', 3, 'color', colors(i,:));
end
legend(legendIds,conditions);
legend boxoff


axis equal
xlim([0 1])

xlabel('restriction')
ylabel('growth')
ylim([0 inf])

%%
figure()
hold all
BCAA = {'HMR_3777' %isoleucine
        'HMR_6923' %leucine
        'HMR_3747'}; %valine
    
labels = {'BCAT1 (isoleucine)', 'BCAT1 (leucine)', 'BCAT1 (valine)'};    
    
for i = 1:length(BCAA)
    curRxn = findIndex(model.rxns, BCAA{i});
    curRxn = ismember(reactionNumbers, curRxn);
    plot(xValues, normalizedGrowth(curRxn,:), 'linewidth', 3, 'color', colors(i,:))
end
plot([0 1], [0,1], 'k--')
legend(labels, 'location', 'se')
legend boxoff


%%
figure()
hold all
conditions = {'glutamate[s]'};
interestingReactions = getBounds(model, conditions);
curRxn = ismember(reactionNumbers, interestingReactions);
plot(xValues, normalizedGrowth(curRxn,:), 'linewidth', 3, 'color', colors(i,:))

plot([0 1], [0,1], 'k--')
legend boxoff