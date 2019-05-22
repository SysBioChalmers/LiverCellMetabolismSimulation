clc
close all
load('model/genericHuman2')
addpath('src')
addpath('displayNetwork')
primColor = [17 115 187]/256;
secColor  = [216 85 39]/256;

epsilon = 10^-6;
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

solution = solveLinMin(model,1);
flux = solution.x * 1000;

[groupNames, reactionGroups, rxnStochiometry, coordinates] = importReactionGroups('rxnGroups/metaboliteMap.txt');

results = zeros(length(groupNames),3);

for i = 1:length(groupNames)
    curRxns = reactionGroups{i};
    curFlux = 0;
    curStoch = rxnStochiometry{i};
    for j = 1:length(curRxns)
        rxnId = findIndex(model.rxns, curRxns{j});
        if not(isempty(rxnId))
            curFlux = curFlux + flux(rxnId) * curStoch(j);
        end
    end
    results(i,1) = curFlux;    
end

if setErrorBounds
    %Set flux error as bounds
    model.ub(findIndex(model.rxns, 'newAlbumin')) = 0;
    fluxError(findIndex(fluxMets,'biomass[s]'),:) = fluxError(findIndex(fluxMets,'biomass[s]'),:)*1000;
    model = bindFBAInterval(model, fluxMets, fluxError/1000);
else
    model.lb(findIndex(model.rxns,'humanGrowhOut')) = -solution.f - epsilon;
    model.ub(findIndex(model.rxns,'humanGrowhOut')) = -solution.f + epsilon;
end

for i = 1:length(reactionGroups)
    curRxns = reactionGroups{i};
    curStoch = rxnStochiometry{i};
    model.c = zeros(length(model.rxns),1);
    
    for j = 1:length(curRxns)
        rxnNr = findIndex(model.rxns, curRxns{j});
        model.c(rxnNr) = curStoch(j);
    end    
    %Maximize
    solution = solveLin(model,1);
    results(i,3) = -solution.f * 1000;
    
    %Minimize
    model.c = -1 * model.c;    
    solution = solveLin(model,1);
    results(i,2) = solution.f * 1000;
end

results(results(:,2)<-1000,2) = -inf;
results(results(:,3)>1000,3) = inf;


%%
networkBg = imread('rxnGroups/network.png');
imshow(networkBg);

for i = 1:length(groupNames)
    fprintf('%s\t%2.2f\n',groupNames{i}, results(i,1))
    XY = coordinates{i};
    if results(i,1)>1
        curStr = sprintf('%1.0f\n(%1.0f, %1.0f)', results(i,1), results(i,2), results(i,3));
    else
        curStr = sprintf('%1.3f\n(%1.3f, %1.3f)', results(i,1), results(i,2), results(i,3));
    end
    text(XY(1),XY(2),curStr, 'Color', secColor, 'FontSize', 10, 'FontName', 'Calibri', 'HorizontalAlignment', 'center')
end
