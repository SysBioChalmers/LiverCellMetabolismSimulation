model = smallModel;

for i = 1:length(model.rxns)
    subRxns = split(model.rxns{i},'+');
    model.rxns{i} = subRxns{1};
end

sources = {'NADPH'};
sinks = {'dGMP', 'dTMP', 'IMP', 'CDP-diacylglycerol', 'cholesterol', 'glycine', 'proline'};
banMets = {'ubiquinol', 'CoA', 'NADP+', 'ubiquinone', 'ATP', 'ADP', 'H+', 'CO2', 'H2O', 'NADH', 'NAD+', 'GMP', 'CTP', 'PPi', 'THF', '5,10-methylene-THF', '5,10-methenyl-THF', '3-phospho-D-glycerate'};
banRxns = {'human_proteinPool'};

[cMatrix, labels, rxnStart, sourceMets, sinkMets] = generateBiPartite(model, smallSolution, sources, sinks, banMets, banRxns);

graph_to_dot('test', cMatrix,  labels, rxnStart, sourceMets, sinkMets) 
