clc
load('model/genericHuman2')
addpath('src')
addpath('displayNetwork')


[fluxMets, fluxValues] = loadFluxes('fluxvalues', 'hepg2-6mm-.txt');
model = setupBiomass(model, 150, 0.5);
model = bindFBA(model, fluxMets, fluxValues(:,2)/1000);


model.lb(findIndex(model.rxns, 'HMR_4940')) = 0;
model.ub(findIndex(model.rxns, 'HMR_6391')) = 0;

model.lb(findIndex(model.rxns, 'HMR_6330')) = 0;
model.ub(findIndex(model.rxns, 'HMR_6330')) = 0;

model.lb(findIndex(model.rxns, 'HMR_6293')) = 0;
model.ub(findIndex(model.rxns, 'HMR_6293')) = 0;

model.lb(findIndex(model.rxns, 'HMR_6289')) = 0;
model.ub(findIndex(model.rxns, 'HMR_6289')) = 0;

model.lb(findIndex(model.rxns, 'HMR_6286')) = 0;
model.ub(findIndex(model.rxns, 'HMR_6286')) = 0;

model.lb(findIndex(model.rxns, 'HMR_4851')) = 0;
model.ub(findIndex(model.rxns, 'HMR_4851')) = 0;


solution = solveLinMin(model,1)


[smallModel, smallSolution] = gemPress(model, solution.x, false, false);
result = balanceValidation(smallModel, smallSolution);


%printExchangeFluxes(model, solution.x);
%printExchangeFluxes(smallModel, smallSolution);

%makeDotGraph(smallModel, smallSolution, {'aspartate[c]', 'asparagine[c]'});
%makeDotGraph(smallModel, smallSolution, {'glutamate[c]', 'L-glutamate 5-semialdehyde[c]'});
%makeDotGraph(smallModel, smallSolution, {'proline[c]', 'L-glutamate 5-semialdehyde[m]'});
%makeDotGraph(smallModel, smallSolution, {'arginine[c]'});
%makeDotGraph(smallModel, smallSolution, {'glucose[c]' , 'glucose-6-phosphate[c]', 'ribose-5-phosphate[c]', 'ribulose-5-phosphate[c]', 'fructose-6-phosphate[c]', 'GAP[c]', 'DHAP[c]', 'D-xylulose-5-phosphate[c]'})

%makeDotGraph(smallModel, smallSolution, {'malate[c]', 'pyruvate[c]', 'L-lactate[c]'})
%makeDotGraph(smallModel, smallSolution, {'5,10-methylene-THF[c]', '5,10-methenyl-THF[c]', 'THF[c]'})
%makeDotGraph(smallModel, smallSolution, {'serine[c]', 'glycine[c]'})
%makeDotGraph(smallModel, smallSolution, {'cysteine[c]','homocysteine[c]', 'methionine[c]'})
%makeDotGraph(smallModel, smallSolution, {'leucine[c]', 'acetoacetyl-CoA[m]', 'cholesterol[c]'});
%makeDotGraph(smallModel, smallSolution, {'(R)-methylmalonyl-CoA[m]', 'propanoyl-CoA[m]', 'pentadecanoyl-[ACP][c]', 'pentadecylic acid[c]', 'margaric acid[c]'})
%makeDotGraph(smallModel, smallSolution, {'(R)-methylmalonyl-CoA[m]', 'propanoyl-CoA[m]', 'succinate[c]'})
%makeDotGraph(smallModel, smallSolution, {'NADPH[c]', 'NADPH[m]'}) 
makeDotGraph(smallModel, smallSolution, {'citrate[c]', 'isocitrate[c]', 'AKG[c]', 'AKG[m]'})

makeDotGraph(smallModel, smallSolution, {'serine[c]', 'glycine[c]', '5,10-methylene-THF[c]', '5,10-methenyl-THF[c]'})
makeDotGraph(smallModel, smallSolution, {'cysteine[c]', 'sulfite[m]', 'sulfite[c]', 'sulfate[c]'})

makeDotGraph(smallModel, smallSolution, {'glutamate[m]', 'glutamate[c]'}) 

%makeDotGraph(smallModel, smallSolution, {'glutamate[m]', 'glutamate[c]', 'aspartate[m]', 'aspartate[c]'}) 

%makeDotGraph(smallModel, smallSolution, {'malate[c]', 'pyruvate[c]', 'L-lactate[c]', 'pyruvate[m]'})
%
%makeDotGraph(smallModel, smallSolution, {'AKG[c]'}) 

%makeDotGraph(smallModel, smallSolution, {'malate[m]','citrate[c]', 'isocitrate[c]', 'OAA[m]','fumarate[m]', 'succinate[m]', 'AKG[m]', 'citrate[c]',  'citrate[m]', 'isocitrate[m]', 'acetyl-CoA[c]', 'acetyl-CoA[m]', 'succinyl-CoA[m]'})
%makeDotGraph(smallModel, smallSolution, {'OAA[m]', 'OAA[c]', 'AKG[m]', 'AKG[c]', 'malate[c]', 'malate[m]'})

%%
%Anaplerosis analysis:
%[smallModel, smallSolution] = gemPress(model, solution.x, true, false);
%[metNames, percents] = printAnaplerosis(smallModel, smallSolution);



