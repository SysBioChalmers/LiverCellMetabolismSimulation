clc
load('model/genericHuman2')
addpath('src')
addpath('displayNetwork')

cellType = 'hepg2';
condition = '22';
fluxPhase = 2;

if and(strcmp(cellType, 'hepg2'),  fluxPhase == 2)
    fileName = ['confidenceIntervalls\output\hepg2-' num2str(condition) '.tsv'];
    raw = IO(fileName);
    fluxMets = raw(2:end,1);
    fluxValues = cell2nummat(raw(2:end,2));
    fluxData = fluxValues/1000;
    fluxError = cell2nummat(raw(2:end,3:4));
%     aroundZero = sign(fluxError(:,1)) == -sign(fluxError(:,2));
%     fluxData(aroundZero) = 0;
else
    [fluxMets, fluxValues] = loadFluxes('fluxvalues', cellType, condition);
    fluxData = fluxValues(:,fluxPhase)/1000;
end

model = setupBiomass(model, 48, 1);

model = bindFBA(model, fluxMets, fluxData);

solution = solveLinMin(model,1);
[smallModel, smallSolution] = gemPress(model, solution.x, true, true);


sources = {'glutamine', 'leucine', 'isoleucine', 'valine', 'phenylalanine', 'arginine', 'cysteine'};
sinks = {'proline', 'glycine', 'asparagine', 'aspartate', 'tyrosine', 'glutamate', 'alanine'};
curencyMets = {'H2O', 'CO2', 'Pi' 'ubiquinol', 'GSH', 'UDP', 'UTP', 'ATP', 'AMP', 'ADP', 'CoA' 'H+',  'NADH', 'NAD+', 'GMP', 'CTP', 'PPi', 'NADP+', 'NADPH', 'THF', '5,10-methylene-THF', '5,10-methenyl-THF'};
poolRxns = {'human_proteinPool', 'metabolitePool'};

%Filter out fluxes with minor contribution to the pool
smallSolution = filterMets(smallModel, smallSolution, [sources sinks], 0.1);

linkAA(smallModel, smallSolution, sources, sinks, curencyMets, poolRxns)





