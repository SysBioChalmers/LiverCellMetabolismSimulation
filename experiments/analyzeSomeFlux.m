clc
load('genericHuman2')
addpath('src')
[id, exchangeRxns] =getExchangeRxns(model);
model = configureSMatrix(model, 40, 'HumanGrowth', 'human_growthMaintainance[c]');
model = configureSMatrix(model, 6, 'HumanGrowth', 'human_protein_pool[c]');
model = configureSMatrix(model, 0.1, 'HumanGrowth', 'human_RNAPool[c]');
model = configureSMatrix(model, 0.1, 'HumanGrowth', 'human_DNAPool[c]');
model = configureSMatrix(model, 1, 'HumanGrowth', 'glycogen[c]');
model = configureSMatrix(model, 0, 'lipidPool', 'fattyAcidPool[c]');
model = configureSMatrix(model, 0.3, 'HumanGrowth', 'lipidPool[c]');


%Set Non-GAM
model = setParam(model, 'lb', 'human_ATPMaintainance', 0.5);
[fluxMets, fluxValues] = loadFluxes('fluxvalues', 'hepg2-0mm-.txt');
nrDists = size(fluxValues,2)-1;

allSolutions = zeros(length(model.rxns), nrDists);
mu = zeros(1,nrDists);
for i = 2:nrDists
    model = bindFBA(model, fluxMets, fluxValues(:,i)/1000);
    objectiveFunction = 'HumanGrowth';
    model = setParam(model, 'obj', objectiveFunction, 1);
    model = setParam(model, 'ub', objectiveFunction, 0.05);
    solution = solveLinMin(model);
    allSolutions(:,i) = solution.x;
    mu(i) = solution.x(findIndex(model.rxns,objectiveFunction));
end
tresh = 10^-6;
results = allSolutions(exchangeRxns,:);

exchange=sum(abs(results),2)>tresh;
results = results(exchange,:);
[res, id] = sort(mean(results,2));
results = results(id,:);

eq = constructEquations(model, exchangeRxns(exchange));
eq = eq(id);
for i = 1:length(eq)
    
   fprintf('%s', eq{i}) 
   for j = 1:length(results(i,:))
       fprintf('\t%f', results(i,j)) 
   end
   fprintf('\n') 
end

%