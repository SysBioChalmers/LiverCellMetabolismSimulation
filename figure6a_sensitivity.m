clc
load('model/genericHuman2')
addpath('src')
addpath('sensitivity analysis')

cellType = 'hepg2';
condition = '22';
setErrorBounds = false;

if strcmp(cellType, 'hepg2')
    fileName = ['confidenceIntervalls\output\hepg2-' num2str(condition) '.tsv'];
    raw = IO(fileName);
    fluxMets = raw(2:end,1);
    fluxValues = cell2nummat(raw(2:end,2));
    fluxData = fluxValues/1000;
    fluxError = cell2nummat(raw(2:end,3:4));
else
    [fluxMets, fluxValues] = loadFluxes('fluxvalues', cellType, condition);
    fluxData = fluxValues(:,2)/1000;
end

model = setupBiomass(model, 48, 1);
model = bindFBA(model, fluxMets, fluxData);

%%%%%%%%%%%%%%%%%%%%
%Relax upper bound for uptake reactions to allow sub maximal flux which is
%required to prevent linear pathways from beeing infeasible
[rxn, id] = getExchangeRxns(model);
for i = 1:length(id)
    if model.lb(id(i)) < 0
        model.ub(id(i)) = 0;
    end
end


curSol = solveLinMin(model,1);
fluxes = curSol.x;
growthTolerance = 0;

tresh = 10^-6;



%%%%%%%%%%%%%%%%%
%Remove reactions that we are not interested in analyzing
toAnalyse = abs(fluxes)>tresh; %only reactions with flux

%ignore exchange rxns
[rxn, id] = getExchangeRxns(model);
toAnalyse(id) = false; 

%Ignore transport reactions and artificial reactions
banedSubs = {'Artificial', 'Putative Reactions', 'Artificall reactions', 'Artificial reactions', 'Transport, mitochondrial', 'Transport, extracellular', 'Transport, nuclear', 'Transport, peroxisomal', 'Transport, lysosomal', 'Transport, lysosome to ER', 'Transport, Golgi apparatus'};
toAnalyse(ismember(model.subSystems, banedSubs)) = false; 

%Constrain objective function 
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


%Sensitivity analysis
confirmedDetermined = find(toAnalyse);
sensitivity = zeros(length(confirmedDetermined),1);

model.c = zeros(length(fluxes),1);
model.c(objRxn) = 1;
model.lb(objRxn) = 0;
model.ub(objRxn) = 1000;

factor = 0.99;
for i = 1:length(confirmedDetermined)
    i
    %100*(1 - (i/nrOfFluxes))
    rxnNr = confirmedDetermined(i);
    
    %save bounds
    lb = model.lb(rxnNr);
    ub = model.ub(rxnNr);
    
    model.lb(rxnNr) = fluxes(rxnNr)*factor;
    model.ub(rxnNr) = fluxes(rxnNr)*factor;
    curSol = solveLin(model,1);
    if length(curSol.x) > 1
        sensitivity(i) = curSol.x(objRxn);
    else
        display('error')
    end
    
    %restore bounds
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

model = mergeSubsystems(model, 'Fatty acid activation (cytosolic)', 'Fatty acid reactions');
model = mergeSubsystems(model, 'Fatty acid biosynthesis (even-chain)', 'Fatty acid reactions');
model = mergeSubsystems(model, 'Fatty acid biosynthesis (odd-chain)', 'Fatty acid reactions');
model = mergeSubsystems(model, 'Fatty acid elongation (even-chain)', 'Fatty acid reactions');
model = mergeSubsystems(model, 'Fatty acid desaturation (even-chain)', 'Fatty acid reactions');
model = mergeSubsystems(model, 'Fatty acid transfer reactions', 'Fatty acid reactions');
model = mergeSubsystems(model, 'Acyl-CoA hydrolysis', 'Fatty acid reactions');

model = mergeSubsystems(model, 'Acyl-CoA hydrolysis', 'Fatty acid reactions');


model = mergeSubsystems(model, 'Amino sugar and nucleotide sugar metabolism', 'Chondroitin / heparan sulfate biosynthesis');


model = mergeSubsystems(model, 'Cholesterol biosynthesis 1 (Bloch pathway) ', 'Cholesterol metabolism');
model = mergeSubsystems(model, 'Terpenoid backbone biosynthesis', 'Cholesterol metabolism');



model = mergeSubsystems(model, 'Purine metabolism', 'Nucleotide metabolism');
model = mergeSubsystems(model, 'Pyrimidine metabolism', 'Nucleotide metabolism');
model = mergeSubsystems(model, 'Tricarboxylic acid cycle and glyoxylate/dicarboxylate metabolism', 'TCA and glyoxylate metabolism');




sensitivityLimit = 0.01;
confirmedSensitive = confirmedDetermined(E>sensitivityLimit);
filteredE = E(E>sensitivityLimit);

determinedSubs = model.subSystems(confirmedSensitive);
subSys = unique(determinedSubs);


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

meanofData = round(meanofData, 4);

T = table(1./meanofData, subSys);
[crap, indx] = sortrows(T, [1 2], 'ascend');

meanofData = meanofData(indx);
data = data(indx);
subSys = subSys(indx);

hold on
color = [93 155 211]/256;

subSysLabel = [subSys];
countSensitive = 0;
countTotal = 0;
tresh = 0.5;

plot([0 length(meanofData)+1], [1 1], 'k-')

plot([0 length(meanofData)+1], tresh * [1 1], 'k-')

for i = 1:length(meanofData)
    bar(i, meanofData(i), 'FaceColor', color,'LineStyle','none');
    nrOfRxns = num2str(length(data{i}));
    subSysLabel{i} = [subSysLabel{i} ' (' nrOfRxns ')'];
    countSensitive = countSensitive + sum(data{i}>tresh);
    countTotal = countTotal + length(data{i});
end

xticks(1:length(meanofData))
xticklabels(subSysLabel)

cordinates = beswarmPlot(data, 0.05);
ylabel('sensitivity')
ylim([0 1.05])


xtickangle(30)

countSensitive
countTotal
