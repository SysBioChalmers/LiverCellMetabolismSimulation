clc
load('model/genericHuman2')
addpath('src')

[fluxMets, fluxValues] = loadFluxes('fluxvalues', 'hepg2-6mm-.txt');
model = setupBiomass(model, 150, 0.5, 0.55);
model = bindFBA(model, fluxMets, fluxValues(:,2)/1000);
solution = solveLinMin(model,1)

tresh = 10^-6;
fluxes = solution.x;
growthTolerance = 0.01;

excludeRxns = abs(fluxes)<tresh;

[rxn, id] = getExchangeRxns(model);
excludeRxns(id) = true; %ignor exchange rxns

objRxn = find(model.c);

%Fix objective function 
model = setParam(model, 'lb', objRxn, fluxes(objRxn)*(1-growthTolerance));
model = setParam(model, 'ub', objRxn, fluxes(objRxn));

%Minimize all fluxes
newObjective = -sign(fluxes);
newObjective(excludeRxns) = 0;
model.c = newObjective;

sum(not(excludeRxns))
curSol = solveLin(model,1);
excludeRxns(abs(curSol.x-fluxes)>10^-6) = true;
sum(not(excludeRxns))

%Maximize all fluxes
newObjective = -sign(fluxes);
newObjective(excludeRxns) = 0;
model.c = newObjective;

curSol = solveLin(model,1);
excludeRxns(abs(curSol.x-fluxes)>10^-6) = true;
sum(not(excludeRxns))

possiblyDeterminedRxns = find(not(excludeRxns));
nrOfFluxes = length(possiblyDeterminedRxns);
for i = nrOfFluxes:-1:1
    %100*(1 - (i/nrOfFluxes))
    variable = false;    
    rxnNr = possiblyDeterminedRxns(i);

    model.c = zeros(length(fluxes),1);
    model.c(rxnNr) = 1;
    curSol = solveLin(model,1);
    if curSol.x(rxnNr) > (fluxes(rxnNr) + tresh)
%         curSol.x(rxnNr)
%         fluxes(rxnNr)
        variable = true;
    else
       model.c(rxnNr) = -1;
       curSol = solveLin(model,1);
       if curSol.x(rxnNr) < (fluxes(rxnNr) - tresh)
%            curSol.x(rxnNr)
%            fluxes(rxnNr)           
           variable = true;          
       end
    end
    
    if variable == true
        possiblyDeterminedRxns(i) = [];
    end
end

%Sensitivity analysis
confirmedDetermined = possiblyDeterminedRxns;
sensitivity = zeros(length(confirmedDetermined),1);

model.c = zeros(length(fluxes),1);
model.c(objRxn) = 1;
model.lb(objRxn) = 0;


factor = 0.99;
for i = 1:length(confirmedDetermined)
    %100*(1 - (i/nrOfFluxes))
    rxnNr = confirmedDetermined(i);
    lb = model.lb(rxnNr);
    ub = model.ub(rxnNr);
    
    model.lb(rxnNr) = fluxes(rxnNr)*factor;
    model.ub(rxnNr) = fluxes(rxnNr)*factor;      
    curSol = solveLin(model,1);
    if length(curSol.x) > 1
        sensitivity(i) = curSol.x(objRxn);
    end
    
    model.lb(rxnNr) = lb;
    model.ub(rxnNr) = ub;    
end
sensitivityNorm = sensitivity/fluxes(objRxn);
E = (1-sensitivityNorm)/(1-factor);

for i = 1:length(confirmedDetermined)
    curRxn = confirmedDetermined(i);
    rxn = constructEquations(model, curRxn);
    sub = model.subSystems{curRxn};
    fprintf('%s\t%s\t%f\t%f\t%s\n', model.rxns{curRxn}, rxn{1}, fluxes(curRxn), E(i), sub)
end

subSystemData = model.subSystems;

% 'Fatty acid biosynthesis (even-chain)'
% 'Fatty acid transfer reactions'
% 'Fatty acid biosynthesis and transfer reactions'
 

determinedSubs = subSystemData(confirmedDetermined);
subSys = unique(determinedSubs);

subSys(ismember(subSys,'Artificial reactions')) = [];
subSys(ismember(subSys,'Artificial')) = [];

% 'Cholesterol biosynthesis 1 (Bloch pathway)'
% 'Cholesterol metabolism'
% 'Cholesterol biosynthesis and metabolism'


fluxSubs = model.subSystems(abs(fluxes)>tresh);
for i = 1:length(subSys)
    detCount = sum(ismember(determinedSubs, subSys{i}));
    fluxCount = sum(ismember(fluxSubs, subSys{i}));
    ratio = detCount/fluxCount;
    fprintf('%s\t%i\t%2.2f\n', subSys{i}, fluxCount, ratio)
end

