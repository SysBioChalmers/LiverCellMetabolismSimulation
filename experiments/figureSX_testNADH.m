clc
load('model/genericHuman2')

results = zeros(3,2);

atpFlux = findIndex(model.rxns, 'human_ATPMaintainance');
cmplx1Flux = findIndex(model.rxns, 'HMR_6921');
relevantReactions = [atpFlux cmplx1Flux];

[id, exchangeRxns] =getExchangeRxns(model);
model = setParam(model, 'lb', exchangeRxns, 0);
model = setParam(model, 'ub', exchangeRxns, 1000);    

upMets = {'H2O[s]', 'O2[s]', 'H+[s]', 'NH3[s]', 'Pi[s]'};
uptake = getBounds(model, upMets);
model = setParam(model, 'lb', uptake, -1000);

% reMets = {'HCO3-[s]', 'CO2[s]', 'H+[s]', 'NH3[s]', 'L-lactate[s]', 'alanine[s]'};
% release = getBounds(model, reMets);
% model = setParam(model, 'ub', release, 1000);

objectiveFunction = 'human_ATPMaintainance';
model = setParam(model, 'obj', relevantReactions, 1 -0.1);
model = setParam(model, 'lb', objectiveFunction, 0);
model = setParam(model, 'ub', objectiveFunction, 1000);


tmpModel = model;
reactionNumbers = getBounds(model, {'L-lactate[s]'});
tmpModel = setParam(tmpModel, 'lb', reactionNumbers, -1);
solution = solveLinMin(tmpModel,1);
results(1,:) = solution.x(relevantReactions);

tmpModel = model;
reactionNumbers = getBounds(tmpModel, {'glutamine[s]'});
tmpModel = setParam(tmpModel, 'lb', reactionNumbers, -1);
solution = solveLinMin(tmpModel,1);
results(2,:) = solution.x(relevantReactions);

tmpModel = model;
reactionNumbers = getBounds(tmpModel, {'glutamine[s]'});
tmpModel = setParam(tmpModel, 'lb', reactionNumbers, -1);
reactionNumbers = getBounds(tmpModel, {'alanine[s]'});
tmpModel = setParam(tmpModel, 'lb', reactionNumbers, 1);
solution = solveLinMin(tmpModel,1);
results(3,:) = solution.x(relevantReactions);

yield = results(:,1)./results(:,2);
labels = {'lactate->CO2', 'glutamine->CO2', 'glutamine->alanine'};

barh(yield)
yticks([1 2 3])
yticklabels(labels)
xlabel('ATP/CI flux [mol/mol]')
xlim([0 4])
