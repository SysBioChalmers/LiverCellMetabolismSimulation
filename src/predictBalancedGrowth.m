function fluxes = predictBalancedGrowth(model, growthRate, mets)
    scaleFactor = 1000;
    metNames = modifyMetNames(model); 
    aminoAcids = model.S(:,findIndex(model.rxns, 'human_proteinPool'));
    biomassRatio = model.S(findIndex(metNames,'human_protein_pool[c]'),findIndex(model.rxns, 'HumanGrowth'));
    aminoAcidId= find(aminoAcids);
    aminoAcidRatio = aminoAcids(aminoAcidId);
    aminoAcidNames = model.metNames(aminoAcidId);
    fluxes = zeros(length(mets),1);
    prefactor = -biomassRatio * growthRate * scaleFactor;
    for i = 1:length(mets)
       namePartA = strsplit(mets{i}, '[');
       namePartA = namePartA{1};
       matchNr = find(ismember(aminoAcidNames, namePartA));
       if ~isempty(matchNr)
            fluxes(i) = aminoAcidRatio(matchNr) * prefactor;
       end
    end
    
end

