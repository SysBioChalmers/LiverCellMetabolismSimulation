clc
load('model/genericHuman2')
addpath('src')
addpath('displayNetwork')


[fluxMets, fluxValues] = loadFluxes('fluxvalues', 'hepg2-0mm-.txt');
model = setupBiomass(model, 150, 0.5);
model = bindFBA(model, fluxMets, fluxValues(:,2)/1000);

solution = solveLinMin(model,1)

[smallModel, smallSolution] = gemPress(model, solution.x, false, false);

%mitochondria
mitoMets = {
'(R)-methylmalonyl-CoA[m]'
'acetoacetyl-CoA[m]'
'acetyl-CoA[m]'
'AKG[m]'
'alanine[m]'
'aspartate[m]'
'citrate[m]'
'fumarate[m]'
'glutamate[m]'
'isocitrate[m]'
'L-glutamate 5-semialdehyde[m]'
'malate[m]'
'OAA[m]'
'proline[c]'
'propanoyl-CoA[m]'
'pyruvate[m]'
'succinate[m]'
'succinyl-CoA[m]'
};

%cytoplasm
cytoMets = {
'3-phosphoserine[c]'
'AKG[c]'
'alanine[c]'
'asparagine[c]'
'aspartate[c]'
'citrate[c]'
'cysteine[c]'
'glutamate[c]'
'glutamine[c]'
'glycine[c]'
'isoleucine[c]'
'leucine[c]'
'lysine[c]'
'malate[c]'
'OAA[c]'
'pyruvate[c]'
'serine[c]'
'valine[c]'
'acetyl-CoA[c]'
'L-lactate[s]'
};

makeMetMetMap(smallModel, smallSolution, mitoMets);

makeMetMetMap(smallModel, smallSolution, cytoMets);

%%
[fluxMets, fluxValues] = loadFluxes('fluxvalues', 'hepg2-6mm-.txt');
model = bindFBA(model, fluxMets, fluxValues(:,3)/1000);
solution = solveLinMin(model,1)

[smallModel, smallSolution] = gemPress(model, solution.x, false, false);

makeMetMetMap(smallModel, smallSolution, mitoMets);

makeMetMetMap(smallModel, smallSolution, cytoMets);
