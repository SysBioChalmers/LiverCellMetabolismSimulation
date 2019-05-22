clc
load('../model/genericHuman2')
addpath('../src')

[id, exchangeRxns] =getExchangeRxns(model);
model = setParam(model, 'lb', exchangeRxns, 0);

model = setParam(model, 'ub', exchangeRxns, 1000);   

upMets = {'H2O[s]', 'O2[s]', 'NH3[s]', 'H+[s]'};
uptake = getBounds(model, upMets);
model = setParam(model, 'lb', uptake, -1000);

% reMets = {'HCO3-[s]', 'CO2[s]', 'H2O[s]', 'L-lactate[s]', 'H+[s]'};
% release = getBounds(model, reMets);
% model = setParam(model, 'ub', release, 1000);

relevantExchange = exchangeRxns;
relevantExchange(ismember(relevantExchange,uptake)) = [];
relevantExchange(ismember(relevantExchange,release)) = [];


%calculate and remove glucose rxns
reactionNumbers = getBounds(model, {'glucose[s]'});
model = setParam(model, 'lb', reactionNumbers, -1);
model = setParam(model, 'ub', reactionNumbers, 0);

results = zeros(length(relevantExchange),1);
for i = 1:length(relevantExchange)
    model = setParam(model, 'obj', relevantExchange(i), 1);
    solution = solveLin(model,1);
    results(i) = -solution.f;
end

glcRxns = relevantExchange(results>10^-6);
relevantExchange(ismember(relevantExchange, glcRxns)) = [];


%calculate cystine results
reactionNumbers = getBounds(model, {'cystine[s]'});
model = setParam(model, 'lb', reactionNumbers, -1);
model = setParam(model, 'ub', reactionNumbers, 0);

results = zeros(length(relevantExchange),1);

for i = 1:length(relevantExchange)
    model = setParam(model, 'obj', relevantExchange(i), 1);
    solution = solveLin(model,1);
    results(i) = -solution.f;
end

cycReactions = relevantExchange(results>10^-6);


fullSolution = zeros(length(cycReactions),length(model.rxns));
for i = 1:length(cycReactions)
    model = setParam(model, 'obj', cycReactions(i), 1);
    solution = solveLinMin(model,1);
    results(i) = -solution.f;
    fullSolution(i,:) = solution.x;
end

%Print all fluxes
for i = 1:size(fullSolution,2)
    if sum(abs(fullSolution(:,i))) > 10^-6
        rxn = constructEquations(model, i);
        fprintf('%s\t%s', model.rxns{i}, rxn{1})
        for j = 1:size(fullSolution,1)
            fprintf('\t%2.2f', fullSolution(j,i))
        end
        fprintf('\n');
    end
end

%List strategies
metList = modifyMetNames(model);
cysMet = ismember(metList, 'cysteine[c]');
cysRxns = find(sum(abs(model.S(cysMet,:))>0,1));
cysRxns(sum(abs(fullSolution(:,cysRxns)),1)<10^-6) = [];



strategies = {'HMR_3912'; %cys 2 mercaptopyruvate/3-mercaptolactate
              'HMR_4302'; %cys 2 H2S 
              'HMR_4324'; %cys 2 GSH
              'HMR_9065'; %cys 2 export
              'cysteineDesulfurase'; %cys 2 S-Sulfanylglutathione
              };    
          
results = zeros(length(cycReactions), length(strategies));

for i = 1:length(strategies)
    tmpModel = model;
    for j = 1:length(strategies)
        if i~=j
            curRxn = findIndex(model.rxns, strategies{j});
            tmpModel.lb(curRxn) = 0;
            tmpModel.ub(curRxn) = 0;
        end
    end
    
    for j = 1:length(cycReactions)
        tmpModel = setParam(tmpModel, 'obj', cycReactions(j), 1);
        solution = solveLin(tmpModel,1);
        results(j, i) = -solution.f;
    end
end

results(abs(results)<10^-6) = 0;


metNames = constructEquations(model, cycReactions);
rxnNames = constructEquations(model, strategies);
tmp = clustergram(results', 'ShowDendrogram', 'off', 'ColumnLabels', metNames, 'RowLabels', rxnNames, 'ColumnLabelsRotate', 90, 'DisplayRatio', [0.1 0.1]);
h = plot(tmp)
h.Position = [0.05 0.55 0.4 0.4];
          