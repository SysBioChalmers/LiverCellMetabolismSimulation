clc
load('model/genericHuman2')
addpath('src')
addpath('sensitivity analysis')

epsilon = 10^-6;
cellType = 'hepg2';
condition = '22';
setErrorBounds = true;
fluxPhase = 2;

constrainedFVA = true;

if and(strcmp(cellType, 'hepg2'),  fluxPhase == 2)
    fileName = ['confidenceIntervalls\output\hepg2-' num2str(condition) '.tsv'];
    raw = IO(fileName);
    fluxMets = raw(2:end,1);
    fluxValues = cell2nummat(raw(2:end,2));
    fluxData = fluxValues/1000;
    fluxError = cell2nummat(raw(2:end,3:4));
else
    [fluxMets, fluxValues] = loadFluxes('fluxvalues', cellType, condition);
    fluxData = fluxValues(:,fluxPhase)/1000;
end

model = setupBiomass(model, 48, 1);

model = bindFBA(model, fluxMets, fluxData);
solution = solveLinMin(model,1)

tresh = 10^-6;
fluxes = solution.x;

%Remove reactions that we are not interested in analyzing
toAnalyse = abs(fluxes)>tresh; %only reactions with flux

if constrainedFVA
    %fix growth rate
    model.lb(model.c==1) = (-solution.f)*(1-epsilon);

    %prevent byproduct formation
    [rxns, id] = getExchangeRxns(model);
    id(abs(fluxes(id))>tresh) = [];
    model.ub(id) = 0;
end

%%%%%%%%%%
%Make structure to store results:
result = zeros(sum(toAnalyse), 4);
result(:,1) = find(toAnalyse); %ids of sensitivity
result(:,2) = fluxes(result(:,1)); % fluxes of sensitivity




%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Run FVA on the uncertain reactions

for i = 1:length(result)
    rxnNr = result(i,1);

    if result(i,3) == 0 
        model.c = zeros(length(fluxes),1);
        model.c(rxnNr) = -1;
        curSol = solveLin(model,1);
        result(i,3) = curSol.x(rxnNr);
    end
    
    if result(i,4) == 0 
        model.c = zeros(length(fluxes),1);
        model.c(rxnNr) = 1;
        curSol = solveLin(model,1);
        result(i,4) = curSol.x(rxnNr);
    end    
end

%%
%%%%%%%%%%%
%Display results

absVar = abs(result(:,4)-result(:,3))/2;
absVarNorm = absVar./abs(result(:,2));
[absVarNorm, indx] = sort(absVarNorm);
result = result(indx,:);



printFVAResult(model, result)

vennResult(result, 0.05)

showSubsystemDetermination(model, result, 0.05)


