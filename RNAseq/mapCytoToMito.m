function result = mapCytoToMito(model)
cyt = find(ismember(model.comps, 'c'));
mit = find(ismember(model.comps, 'm'));

allMitRxns = model.rxnComps == mit;
allMitRxns(contains(model.subSystems, 'Transport')) = 0;
allMitRxns = find(allMitRxns);

%same rxn if same metabolites involved
binaryS = sign(abs(model.S));
cytosolicMets = model.metComps == cyt;

mappedRxns = zeros(size(binaryS,1), length(allMitRxns));

for i = 1:length(allMitRxns)
    currentRxn = binaryS(:,allMitRxns(i));
    currentMets = find(currentRxn);
    
    %map metabolites
    for j = 1:length(currentMets)
        metName = model.metNames{currentMets(j)};
        mapedMet = and(ismember(model.metNames,metName), cytosolicMets);
        mappedRxns(mapedMet,i) = 1;
    end
   
end

[~,index_A,index_B] = intersect(mappedRxns',binaryS','rows');

mitRxns = allMitRxns(index_A);
cytRxns = index_B;
result = [cytRxns mitRxns];
end
