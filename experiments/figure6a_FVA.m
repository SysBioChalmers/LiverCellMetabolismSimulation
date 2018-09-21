clc
load('../model/genericHuman2')
addpath('../src')
addpath('../sensitivity analysis')



[fluxMets, fluxValues] = loadFluxes('../fluxvalues', 'hepg2-6mm-.txt');
model = setupBiomass(model, 150, 0.5);

model = bindFBA(model, fluxMets, fluxValues(:,2)/1000);
solution = solveLinMin(model,1)

tresh = 10^-6;
fluxes = solution.x;
growthTolerance = 0.01;



%Remove reactions that we are not interested in analyzing
toAnalyse = abs(fluxes)>tresh; %only reactions with flux


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


