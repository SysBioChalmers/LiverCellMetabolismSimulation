load('../model/genericHuman2')
tresh = 10^-6;

model = configureSMatrix(model, 160, 'HumanGrowth', 'human_growthMaintainance[c]');
model = configureSMatrix(model, 6, 'HumanGrowth', 'human_protein_pool[c]');
model = configureSMatrix(model, 0.1, 'HumanGrowth', 'human_RNAPool[c]');
model = configureSMatrix(model, 0.1, 'HumanGrowth', 'human_DNAPool[c]');
model = configureSMatrix(model, 1, 'HumanGrowth', 'glycogen[c]');
model = configureSMatrix(model, 0, 'lipidPool', 'fattyAcidPool[c]');
model = configureSMatrix(model, 0.3, 'HumanGrowth', 'lipidPool[c]');
model = configureSMatrix(model, 0.1118 	* 0.047/(0.057+0.047), 'human_proteinPool', 'glutamine[c]');
model = configureSMatrix(model, 0.1118 	* 0.057/(0.057+0.047), 'human_proteinPool', 'glutamate[c]');

model = setParam(model, 'lb', 'human_ATPMaintainance', 0.5);

objectiveFunction = 'HumanGrowth';
model = setParam(model, 'obj', objectiveFunction, 1);
model = setParam(model, 'lb', objectiveFunction, 0);
model = setParam(model, 'ub', objectiveFunction, 1000);

%ExchangeFluxes A
[fluxMetsA, fluxValuesA] = loadFluxes('../fluxvalues', 'hepg2-6mm-.txt');
modelA = bindFBA(model, fluxMetsA, fluxValuesA(:,2)/1000);
pFBAsolutionA = solveLinMin(modelA,1);
pFBAA = pFBAsolutionA.x;

%ExchangeFluxes B
[fluxMetsB, fluxValuesB] = loadFluxes('../fluxvalues', 'hepg2-0mm-.txt');
modelB = bindFBA(model, fluxMetsB, fluxValuesB(:,3)/1000);
pFBAsolutionB = solveLinMin(modelB,1);
pFBAB = pFBAsolutionB.x;

%constrain growth
modelA = setParam(modelA, 'lb', objectiveFunction, -pFBAsolutionA.f);
modelA = setParam(modelA, 'ub', objectiveFunction, -pFBAsolutionA.f);

modelB = setParam(modelB, 'lb', objectiveFunction, -pFBAsolutionB.f);
modelB = setParam(modelB, 'ub', objectiveFunction, -pFBAsolutionB.f);


allFlux = [pFBAA pFBAB];

for i = 1:length(allFlux)   
    if xor(abs(allFlux(i,1))>tresh, abs(allFlux(i,2))>tresh)
        rxn = model.rxns{i};
        eq = constructEquations(model, i);
        sub = model.subSystems{i};
        flux = 1000*allFlux(i,:);
        fprintf('%s\t%s\t%s\t%2.2f\t%2.2f\n', rxn, eq{1}, sub, flux(1), flux(2));
    end
end

%constrain consensus reactions
consensus = find(and(abs(allFlux(:,1))>tresh, abs(allFlux(:,2))>tresh));

modelA = setParam(modelA, 'lb', consensus, allFlux(consensus,1));
modelA = setParam(modelA, 'ub', consensus, allFlux(consensus,1));

modelB = setParam(modelB, 'lb', consensus, allFlux(consensus,2));
modelB = setParam(modelB, 'ub', consensus, allFlux(consensus,2));


%run qMOMA
%Whilst not allowing crazy values 
[momaA, momaB, flag]=qMOMA(modelA,modelB, 2);

allFlux = [pFBAA momaA pFBAB momaB];

fprintf('\n\n------------\n')

%Print differences
for i = 1:length(allFlux)
    if xor(abs(allFlux(i,2))>tresh, abs(allFlux(i,4))>tresh)
        rxn = model.rxns{i};
        eq = constructEquations(model, i);
        sub = model.subSystems{i};
        flux = 1000*allFlux(i,:);
        fprintf('%s\t%s\t%s\t%2.2f\t%2.2f\t%2.2f\t%2.2f\n', rxn, eq{1}, sub, flux(1), flux(2), flux(3), flux(4));
    end
end




