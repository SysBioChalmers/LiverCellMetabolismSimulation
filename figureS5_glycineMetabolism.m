clc
load('model/genericHuman2')
addpath('src')
addpath('displayNetwork')

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

% includedMets = {
%     'glycine[c]'
%     'serine[c]'
%     '3-phospho-D-glycerate[c]'
%     '3-phosphoserine[c]'
%     '3-phosphonooxypyruvate[c]'
%     '3-phosphoserine[c]'
%     };
% makeMetMetMap(model, solution.x, includedMets);

includedMets = {
    'glycine[c]'
    'glycine[m]'
    'D-xylulose-5-phosphate[c]'
    'DHAP[c]'
    'glyoxalate[m]'
    'glyoxalate[c]'
    'glycolate[c]'
    'glycolaldehyde[c]'
    'D-xylulose-1-phosphate[c]'
    'D-xylulose[c]'
    'DHAP[c]'
    };
involvedRxns = makeMetMetMap(model, solution.x, includedMets);
eqns = constructEquations(model, involvedRxns)
genes = model.grRules(involvedRxns)


%[smallModel, smallSolution] = gemPress(model, solution.x, false, false);
%result = balanceValidation(smallModel, smallSolution);
%makeDotGraph(smallModel, smallSolution, {'glycine[c]'}) 

