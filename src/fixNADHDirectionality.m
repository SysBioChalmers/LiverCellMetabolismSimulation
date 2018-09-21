function model = fixNADHDirectionality(model)
NADPHReactions = {
    'HMR_4586'
    'HMR_3804'
    'HMR_4687'
    'HMR_4442'
    'HMR_4333'
    'HMR_4345'
    'HMR_8709'
    'HMR_4590'
    'HMR_8698'
    'HMR_4655'
    'HMR_3995'
    'HMR_3855'
    'HMR_8509'
    'HMR_4686'
    'HMR_8699'
    'HMR_1329'
    'HMR_4687'
    'HMR_4785'
    'HMR_4523'
    'HMR_8344'
    'HMR_6537'
    'HMR_7942'
    'HMR_4534'
    };

NADHReactions = {
    'HMR_4588'
    'HMR_3802'
    'HMR_8144'
    'HMR_4332'
    'HMR_4344'
    'HMR_3925'
    'HMR_6638'
    'HMR_4591'
    'HMR_6646'
    'HMR_4654'
    'HMR_3996'
    'HMR_8508'
    'HMR_8507'
    'HMR_4683'
    'HMR_6675'
    'HMR_1325'
    'HMR_4685'
    'HMR_8097'
    'HMR_4448'
    'HMR_7943'
    'HMR_4531'
    'HMR_3869'
};

model = addConstraint(model, NADPHReactions, 'NADPH');
model = addConstraint(model, NADHReactions, 'NAD+');
end

function model = addConstraint(model, rxns, definingSubstrate)
    metaboliteMap = ismember(model.metNames, definingSubstrate);
    for i = 1:length(rxns)
        rxnNr = findIndex(model.rxns, rxns{i});
        %constructEquations(model, rxnNr)
        direction = full(model.S(metaboliteMap,rxnNr));
        direction = direction(find(direction));
        if direction == -1
            model.lb(rxnNr) = 0;
        else
            model.ub(rxnNr) = 0;
        end
    end
end
