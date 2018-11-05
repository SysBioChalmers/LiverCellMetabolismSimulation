clc
load('model/genericHuman2')
addpath('src')
addpath('displayNetwork')

[fluxMets, fluxValues] = loadFluxes('fluxvalues', 'hepg2-6mm-.txt');

model = setupBiomass(model, 150, 0.5);

model = bindFBA(model, fluxMets, fluxValues(:,2)/1000);

solution = solveLinMin(model,1);
[smallModel, smallSolution] = gemPress(model, solution.x, true, true);


sources = {'glutamine', 'leucine', 'isoleucine', 'valine', 'phenylalanine', 'arginine', 'cysteine'};
sinks = {'proline', 'glycine', 'asparagine', 'aspartate', 'tyrosine', 'glutamate', 'alanine', 'serine'};
curencyMets = {'H2O', 'CO2', 'Pi' 'ubiquinol', 'GSH', 'UDP', 'UTP', 'ATP', 'AMP', 'ADP', 'CoA' 'H+',  'NADH', 'NAD+', 'GMP', 'CTP', 'PPi', 'NADP+', 'NADPH', 'THF', '5,10-methylene-THF', '5,10-methenyl-THF'};
poolRxns = {'human_proteinPool', 'metabolitePool'};

%Filter out fluxes with insignificant contribution to the pool
smallSolution = filterMets(smallModel, smallSolution, [sources sinks], 0.15);


linkAA(smallModel, smallSolution, sources, sinks, curencyMets, poolRxns)





