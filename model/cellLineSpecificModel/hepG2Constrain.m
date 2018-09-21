function model = hepG2Constrain(model)
    kockRxns = importdata('hepKnock.txt');

    affectedRxns = ismember(model.rxns, kockRxns);
    model.lb(affectedRxns) = 0;
    model.ub(affectedRxns) = 0;
    
end