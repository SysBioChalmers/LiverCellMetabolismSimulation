function displayReactionComponents(model, rxn)
    objectiveFunction = findIndex(model.rxns, rxn);
    bioMets = find(model.S(:,objectiveFunction));
    metNames = modifyMetNames(model);
    stochiometry = full(model.S(bioMets,objectiveFunction));

    sumOfMass = 0;
    for i = 1:length(bioMets)
        name = metNames{bioMets(i)};
        s = -stochiometry(i)/1000; %mmol
        molWeight = getMolecularWeight(model, name);
        mass = s*molWeight;
        
        if mass > 0
            sumOfMass = sumOfMass + mass;
            fprintf('%s\t%2.5f\t%2.5f\t%2.5f\n',name, 1000*s, molWeight, mass);
        end
    end
    fprintf('total\t\t\t%2.5f\n', sumOfMass);
end