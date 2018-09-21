clc

tresh = 10^-6;
tempModel = model;
fluxes = solution.x;
growthTolerance = 0.01;


variability = abs(fluxes)<tresh;

[rxn, id] = getExchangeRxns(tempModel);
variability(id) = true; %ignor exchange rxns

objRxn = find(tempModel.c);

%Fix objective function 
tempModel = setParam(tempModel, 'lb', objRxn, fluxes(objRxn)*(1-growthTolerance));
tempModel = setParam(tempModel, 'ub', objRxn, fluxes(objRxn));

%Minimize all fluxes
newObjective = -sign(fluxes);
newObjective(variability) = 0;
tempModel.c = newObjective;

sum(not(variability))
curSol = solveLin(tempModel,1);
variability(abs(curSol.x-fluxes)>10^-6) = true;
sum(not(variability))

%Maximize all fluxes
newObjective = -sign(fluxes);
newObjective(variability) = 0;
tempModel.c = newObjective;

curSol = solveLin(tempModel,1);
variability(abs(curSol.x-fluxes)>10^-6) = true;
sum(not(variability))

possiblyDeterminedRxns = find(not(variability));

%Sensitivity analysis
confirmedDetermined = possiblyDeterminedRxns;
sensitivity = zeros(length(confirmedDetermined),1);

tempModel.c = zeros(length(fluxes),1);
tempModel.c(objRxn) = 1;
tempModel.lb(objRxn) = 0;


factor = 0.99;
for i = 1:length(confirmedDetermined)
    %100*(1 - (i/nrOfFluxes))
    rxnNr = confirmedDetermined(i);
    lb = tempModel.lb(rxnNr);
    ub = tempModel.ub(rxnNr);
    
    tempModel.lb(rxnNr) = fluxes(rxnNr)*factor;
    tempModel.ub(rxnNr) = fluxes(rxnNr)*factor;      
    curSol = solveLin(tempModel,1);
    if length(curSol.x) > 1
        sensitivity(i) = curSol.x(objRxn);
    end
    
    tempModel.lb(rxnNr) = lb;
    tempModel.ub(rxnNr) = ub;    
end
sensitivityNorm = sensitivity/fluxes(objRxn);
E = (1-sensitivityNorm)/(1-factor);

for i = 1:length(confirmedDetermined)
    curRxn = confirmedDetermined(i);
    rxn = constructEquations(tempModel, curRxn);
    sub = tempModel.subSystems{curRxn};
    fprintf('%s\t%s\t%f\t%f\t%s\n', tempModel.rxns{curRxn}, rxn{1}, fluxes(curRxn), E(i), sub)
end

subSystemData = tempModel.subSystems;
subSystemData = mergeSubsystems(subSystemData, 'Fatty acid biosynthesis (even-chain)', 'Fatty acid transfer reactions', 'Fatty acid biosynthesis and transfer reactions');
subSystemData = mergeSubsystems(subSystemData, 'Cholesterol biosynthesis 1 (Bloch pathway) ', 'Cholesterol metabolism', 'Cholesterol biosynthesis and metabolism');
subSystemData = mergeSubsystems(subSystemData, 'Purine metabolism', 'Pyrimidine metabolism', 'Nucleotide metabolism');


sensitivityLimit = 0.01;
confirmedSensitive = confirmedDetermined(E>sensitivityLimit);
filteredE = E(E>sensitivityLimit);

determinedSubs = subSystemData(confirmedSensitive);
subSys = unique(determinedSubs);

subSys(ismember(subSys,'Artificial reactions')) = [];
subSys(ismember(subSys,'Artificial')) = [];



for i = 1:length(subSys)
    affectedRxns = ismember(determinedSubs, subSys{i});
    rxnCount = sum(affectedRxns);
    rxnIds = confirmedSensitive(affectedRxns);
    geneData = tempModel.rxnGeneMat(rxnIds,:);
    uniqueRxnGeneCombo = unique(geneData,'rows');
    %affectedGenes = tempModel.genes(geneData>0);
    %geneCount = length(affectedGenes);
    geneCount = size(uniqueRxnGeneCombo,1);
    avgE = mean(filteredE(ismember(determinedSubs, subSys{i})));
    fprintf('%s\t%i\t%i\t%2.2f\n', subSys{i}, geneCount, rxnCount, avgE)
end

%%
data = cell(length(subSys),1);
meanofData = zeros(length(subSys),1);

for i = 1:length(subSys)
    affectedRxns = ismember(determinedSubs, subSys{i});
    affectedIndx = find(affectedRxns);
    rxnIds = confirmedSensitive(affectedRxns);
    geneData = tempModel.rxnGeneMat(rxnIds,:);
    [uniqueRxnGeneCombo, indx] = unique(geneData,'rows');
    data{i} = filteredE(affectedIndx(indx));
    meanofData(i) = mean(data{i});
end

[meanofData, indx] = sort(meanofData, 1);
data = data(indx);
subSys = subSys(indx);

hold on
color = [93 155 211]/256;
for i = 1:length(meanofData)
    barh(i, meanofData(i), 'FaceColor', color,'LineStyle','none');
    barText = sprintf('%2.2f', avgE);
end
subSysLabel = [subSys];
yticks(1:length(meanofData))
yticklabels(subSysLabel)

cordinates = beswarmPlot(data, 0.08);
