clc
load('model/genericHuman2')
addpath('src')
addpath('displayNetwork')

cellType = 'hepg2';
condition = '22';

if strcmp(cellType, 'hepg2')
    fileName = ['confidenceIntervalls\output\hepg2-' num2str(condition) '.tsv'];
    raw = IO(fileName);
    fluxMets = raw(2:end,1);
    fluxValues = cell2nummat(raw(2:end,2));
    fluxData = fluxValues/1000;
    fluxError = cell2nummat(raw(2:end,3));
else
    [fluxMets, fluxValues] = loadFluxes('fluxvalues', cellType, condition);
    fluxData = fluxValues(:,2)/1000;
end

model = setupBiomass(model, 50, 0.5);
model = bindFBA(model, fluxMets, fluxData);

solution = solveLinMin(model,1);

[smallModel, smallSolution] = gemPress(model, solution.x, false, false);

result = balanceValidation(smallModel, smallSolution);

includedMets = {
'(R)-methylmalonyl-CoA[m]'
'3-phosphoserine[c]'
'acetoacetyl-CoA[m]'
'acetyl-CoA[m]'
'AKG[c]'
'AKG[m]'
'alanine[c]'
'alanine[m]'
'asparagine[c]'
'aspartate[c]'
'aspartate[m]'
'citrate[c]'
'citrate[m]'
'cysteine[c]'
'fumarate[m]'
'glutamate[c]'
'glutamate[m]'
'glutamine[c]'
'glycine[c]'
'isocitrate[m]'
'isoleucine[c]'
'leucine[c]'
'L-glutamate 5-semialdehyde[m]'
'malate[c]'
'malate[m]'
'OAA[c]'
'OAA[m]'
'proline[c]'
'propanoyl-CoA[m]'
'pyruvate[c]'
'serine[c]'
'succinate[m]'
'succinyl-CoA[m]'
'valine[c]'};

makeMetMetMap(smallModel, smallSolution, includedMets);


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

%
includedMets = {'phosphopantetheine[m]', 'dephospho-CoA[m]', 'pantetheine[m]', 'pantothenate[c]', 'pantetheine[c]', 'cysteamine[c]', 'dephospho-CoA[m]', 'CoA[m]', 'hypotaurine[c]', 'taurine[c]'};
makeMetMetMap(model, 1000*solution.x, includedMets);

includedMets = {'(R)-4-phosphopantothenoyl-cysteine[c]', 'D-4-phosphopantothenate[c]', 'phosphopantetheine[c]', 'dephospho-CoA[c]','pantothenate[c]', 'CoA[c]', 'cysteine[c]'};
makeMetMetMap(model, 1000*solution.x, includedMets);
