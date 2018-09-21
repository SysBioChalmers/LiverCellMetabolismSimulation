model = smallModel;

for i = 1:length(model.rxns)
    subRxns = split(model.rxns{i},'+');
    model.rxns{i} = subRxns{1};
end

sources = {'glucose', 'glutamine'};
sinks = {'glutamate', 'alanine'};
banMets = {'ATP', 'ADP', 'H+', 'H2O', 'NADH', 'NAD+', 'GMP', 'CTP', 'PPi', 'NADP+', 'NADPH', 'THF', '5,10-methylene-THF', '5,10-methenyl-THF', '3-phospho-D-glycerate'};
banRxns = {'HMR_4396'};

[cMatrix, labels, rxnStart, sourceMets, sinkMets] = generateBiPartite(model, smallSolution, sources, sinks, banMets, banRxns);

graph_to_dot('test', cMatrix,  labels, rxnStart, sourceMets, sinkMets) 
