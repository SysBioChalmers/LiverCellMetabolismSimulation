clc
load('../model/genericHuman2')
addpath('../src')

[fluxMets, fluxValues] = loadFluxes('../fluxvalues', 'hepg2-6mm-.txt');

model = setupBiomass(model, 150, 0.5, 0.83);

%Set Non-GAM

model = bindFBA(model, fluxMets, fluxValues(:,2)/1000);
reactionNumbers = getBounds(model, {'glutamine[s]', 'glutamate[s]', 'alanine[s]'});

growthRates = zeros(4,1);

%Reference
solution = solveLin(model,1);
growthRates(1) = -solution.f;

%No glutamine
model = setParam(model, 'lb', reactionNumbers(1), 0);
model = setParam(model, 'ub', reactionNumbers(1), 0);
solution = solveLin(model,1);
growthRates(2) = -solution.f;

%No glutamine or glutamate production
model = setParam(model, 'lb', reactionNumbers(2), 0);
model = setParam(model, 'ub', reactionNumbers(2), 0);
solution = solveLin(model,1);
growthRates(3) = -solution.f;

%No glutamine or glutamate production
model = setParam(model, 'lb', reactionNumbers(3), 0);
model = setParam(model, 'ub', reactionNumbers(3), 0);
solution = solveLin(model,1);
growthRates(4) = -solution.f;

growthRates

growthRates./growthRates(1)
%%
load('../model/genericHuman2')
model = setupBiomass(model, 150, 0.5, 0);

%Set Non-GAM

model = bindFBA(model, fluxMets, fluxValues(:,2)/1000);
reactionNumbers = getBounds(model, {'glucose[s]', 'L-lactate[s]'});


growthRates = zeros(3,1);

%Reference
solution = solveLin(model,1);
growthRates(1) = -solution.f;

%No glucose
model = setParam(model, 'lb', reactionNumbers(1), 0);
model = setParam(model, 'ub', reactionNumbers(1), 0);
solution = solveLin(model,1);
growthRates(2) = -solution.f;

%No glucose or lactate production
model = setParam(model, 'lb', reactionNumbers(2), 0);
model = setParam(model, 'ub', reactionNumbers(2), 0);
solution = solveLin(model,1);
growthRates(3) = -solution.f;

growthRates

growthRates./growthRates(1)
%%
load('../model/genericHuman2')
model = setupBiomass(model, 150, 0.5, 0.83);
model = bindFBA(model, fluxMets, fluxValues(:,2)/1000);

growthRates = zeros(4,1);

%Reference
solution = solveLin(model,1);
growthRates(1) = -solution.f;

%No Complex 1
model = setParam(model, 'ub', 'HMR_6921', 0);
solution = solveLin(model,1);
growthRates(2) = -solution.f;

%No ATP synthase
model = setParam(model, 'ub', 'HMR_6916', 0);
solution = solveLinMin(model,1);
growthRates(3) = -solution.f;

%No COX
model = setParam(model, 'ub', 'HMR_6914', 0);
solution = solveLin(model,1);
growthRates(4) = -solution.f;
growthRates

growthRates./growthRates(1)


%%
load('../model/genericHuman2')
model = setupBiomass(model, 150, 0.5, 0.83);
model = bindFBA(model, fluxMets, fluxValues(:,2)/1000);

growthRates = zeros(2,1);

%Reference
solution = solveLin(model,1);
growthRates(1) = -solution.f;

%Allow urea
model.ub(findIndex(model.rxns, 'HMR_4949')) = 1000;
model.lb(findIndex(model.rxns, 'HMR_3809')) = -1000;
model.ub(findIndex(model.rxns, 'HMR_3809')) = 1000;


%Force Urea
model = setParam(model, 'ub', 'HMR_9073', 0); %NH3[x]

[crap, exId] = getExchangeRxns(model);

solution = solveLin(model,1);
growthRates(2) = -solution.f;

growthRates
growthRates./growthRates(1)

exchangeFlux = solution.x(exId);
exId = exId(exchangeFlux>10^-6);
exId = exId(model.ub(exId)>500);

tmp = solution.x(exId);
solution = zeros(length(solution.x),1);
solution(exId) = tmp;
printExchangeFluxes(model, solution)

%%
load('../model/genericHuman2')
model = setupBiomass(model, 150, 0.5, 0.83);
model = bindFBA(model, fluxMets, fluxValues(:,2)/1000);

growthRates = zeros(3,1);

%Reference
solution = solveLin(model,1);
growthRates(1) = -solution.f;

%force palmitate production
model = setParam(model, 'lb', 'EXC_BOTH_m02674s', 0.66*6.3*10^-3); %NH3[x]

[crap, exId] = getExchangeRxns(model);

solution = solveLin(model,1);
growthRates(2) = -solution.f;


model = setParam(model, 'lb', 'EXC_BOTH_m02674s', 0.66*112*10^-3); %NH3[x]

[crap, exId] = getExchangeRxns(model);

solution = solveLin(model,1);
growthRates(3) = -solution.f;


growthRates
growthRates./growthRates(1)


