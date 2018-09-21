clc
tresh = 10^-6;
[rxn, id] = getExchangeRxns(model);
for i = 1:length(id)
    if model.lb(id(i)) < 0
        model.ub(id(i)) = max(0, model.ub(id(i)));
    end
    if model.ub(id(i)) < 1000
        model.ub(id(i)) = 1000;
    end
end

%Constrain other fluxes
% model.ub(solution.x>tresh) = solution.x(solution.x>tresh);
% model.lb(solution.x<-tresh) = solution.x(solution.x<-tresh);


curSol = solveLin(model,1);



fluxes = curSol.x;
growthTolerance = 0;




variability = abs(fluxes)<tresh;

[rxn, id] = getExchangeRxns(model);
variability(id) = true; %ignor exchange rxns

objRxn = find(model.c);

%Fix objective function 
model = setParam(model, 'lb', objRxn, fluxes(objRxn)*(1-growthTolerance));
model = setParam(model, 'ub', objRxn, fluxes(objRxn));

%Minimize all fluxes
newObjective = -sign(fluxes);
newObjective(variability) = 0;
model.c = newObjective;

sum(not(variability))
curSol = solveLin(model,1);
variability(abs(curSol.x-fluxes)>10^-6) = true;
sum(not(variability))

%Maximize all fluxes
newObjective = -sign(fluxes);
newObjective(variability) = 0;
model.c = newObjective;

curSol = solveLin(model,1);
variability(abs(curSol.x-fluxes)>10^-6) = true;
sum(not(variability))

possiblyDeterminedRxns = find(not(variability));

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
subSystemData = mergeSubsystems(subSystemData, 'Fatty acid biosynthesis (even-chain)', 'Fatty acid transfer reactions', 'Fatty acid biosynthesis and transfer reactions');
subSystemData = mergeSubsystems(subSystemData, 'Cholesterol biosynthesis 1 (Bloch pathway) ', 'Cholesterol metabolism', 'Cholesterol biosynthesis and metabolism');
subSystemData = mergeSubsystems(subSystemData, 'Purine metabolism', 'Pyrimidine metabolism', 'Nucleotide metabolism');
subSystemData(ismember(subSystemData,'Tricarboxylic acid cycle and glyoxylate/dicarboxylate metabolism')) = {'TCA and glyoxylate metabolism'};


sensitivityLimit = 0.01;
confirmedSensitive = confirmedDetermined(E>sensitivityLimit);
filteredE = E(E>sensitivityLimit);



determinedSubs = subSystemData(confirmedSensitive);
subSys = unique(determinedSubs);

subSys(ismember(subSys,'Artificial reactions')) = [];
subSys(ismember(subSys,'Artificial')) = [];
subSys(ismember(subSys,'Transport, mitochondrial')) = [];
subSys(ismember(subSys,'Transport, extracellular')) = [];


%%
data = cell(length(subSys),1);
meanofData = zeros(length(subSys),1);

filteredE(filteredE>1) = 1;

for i = 1:length(subSys)
    affectedRxns = ismember(determinedSubs, subSys{i});
    affectedIndx = find(affectedRxns);
    rxnIds = confirmedSensitive(affectedRxns);
    geneData = model.rxnGeneMat(rxnIds,:);
    [uniqueRxnGeneCombo, indx] = unique(geneData,'rows');
    data{i} = filteredE(affectedIndx(indx));
    meanofData(i) = mean(data{i});
end

[meanofData, indx] = sort(meanofData, 1);
data = data(indx);
subSys = subSys(indx);

hold on
color = [93 155 211]/256;

subSysLabel = [subSys];
for i = 1:length(meanofData)
    barh(i, meanofData(i), 'FaceColor', color,'LineStyle','none');
    nrOfRxns = num2str(length(data{i}));
    subSysLabel{i} = [subSysLabel{i} ' (' nrOfRxns ')'];
end

yticks(1:length(meanofData))
yticklabels(subSysLabel)

cordinates = beswarmPlot(data, 0.05);
xlabel('sensitivity')
xlim([0 1.05])
