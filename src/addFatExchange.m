function model = addFatExchange(model, compartment, subSystemName)
% addFatExchange
%   Adds exchange reactions for the metabolites hardcoded in this function
%   using RAVENS frunctionality "addExchangeRxns"
%
%   model            a model structure
%   compartment      the string abbreviation of a compartment e.g. 's'
%   subSystemName    a string with the name of the subsystem e.g. 'fat
%                    rxns'
%   model            an updated model structure
%
%   Avlant Nilsson, 2016-05-16
%


     labels = {
    'decanoic acid'
    'lauric acid'
    'myristic acid'
    'pentadecylic acid'
    'margaric acid'
    'palmitate'
    'stearate'
    'eicosanoate'
    'behenic acid'
    'lignocerate'
    'physeteric acid'
    'palmitolate'
    'oleate'
    'cis-vaccenic acid'
    'cis-gondoic acid'
    'cis-erucic acid'
    'gamma-linolenate'
    '(11Z,14Z)-eicosadienoic acid'
    'dihomo-gamma-linolenate'
    'arachidonate'
    'adrenic acid'
    'EPA'
    'DPA'
    'DHA'
    };

    labels = strcat(labels, ['[' compartment ']']);

    modMetName = modifyMetNames(model);
    metaboliteNumbers = getIndexFromText(modMetName, labels);

    %Test if metabolites exist in model, output errors
    if sum(metaboliteNumbers==-1) ~= 0
        disp('problem with:');
        problemList = labels(metaboliteNumbers==-1)
        metaboliteNumbers(metaboliteNumbers==-1) = [];
    end
    
    %Remove reactions that allready exist
    [exchangeRxns,exchangeRxnsIndexes] = getExchangeRxns(model,'both');
    existingMetId = find(sum(model.S(:,exchangeRxnsIndexes),2));    
    metaboliteNumbers = setdiff(metaboliteNumbers, existingMetId);
    mets = model.mets(metaboliteNumbers);
    
    %Add the exchange reactions
    [model, addedRxns] = addExchangeRxns(model,'both',mets);
    
    %Modify the subsystem name
    addedRxns = find(ismember(model.rxns, addedRxns));
    for i = 1:length(addedRxns)
        model.subSystems{addedRxns(i)} = subSystemName;
    end
end

