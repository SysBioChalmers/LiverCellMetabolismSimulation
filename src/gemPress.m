function [smallModel, smallSolution] = gemPress(model, fullSolution, rmComps, rmDuplicates)
% gemPress
% Simplifies the model by joining linear pathways to lumped reactions
%
%   model               a model struct
%   fullSolution        a vector with fluxes
%   rmComps             boolean, remove compartments?
%   rmDuplicates        boolean, remove duplicated reactions?
%
%   smallModel          reduced model
%   smallSolution       reduced flux vector
%
%   Avlant Nilsson, 2016-05-17
%
[smallModel, smallSolution] = simplifyFluxResult(model, fullSolution, rmComps, rmDuplicates);

nrMets = size(smallModel.S, 1);
nrRxns = size(smallModel.S, 2);

while true
    [smallModel, smallSolution] = simplifyFluxResult(smallModel, smallSolution, false, false);
    
    if and(size(smallModel.S, 1)==nrMets, size(smallModel.S, 2)==nrRxns)
       break 
    else
        nrMets = size(smallModel.S, 1);
        nrRxns = size(smallModel.S, 2);
    end
end

end

function [newModel, fullSolution] = simplifyFluxResult(model, fullSolution, rmComps, rmDuplicates)

fluxThresh = 10^-7;

if rmComps
    [newModel, fullSolution] = removeCompartments(model, fullSolution);
else
    newModel = model;
end

removeRxns = abs(fullSolution)<fluxThresh;

newModel = removeIndicatedRxns(newModel, removeRxns);
fullSolution(removeRxns) = [];

removeMets = sum(abs(newModel.S),2)==0;
newModel = removeIndicatedMets(newModel, removeMets);


%Make model positive
[newModel, fullSolution] = makePositiveModel(newModel, fullSolution);

[id, exchangeRxns] = getExchangeRxns(newModel);
exchangeMets = find(sum(abs(newModel.S(:,exchangeRxns)),2));

%simplify all metabolites with a reaction Participation = 2
    for i=length(newModel.mets):-1:1
        if not(ismember(i, exchangeMets))
           rxnParticipation = find(newModel.S(i,:));
           if length(rxnParticipation) == 2
               stochA = newModel.S(i, rxnParticipation(1));
               stochB = newModel.S(i, rxnParticipation(2));

               if areIntegers(stochA, stochB)

                  if abs(stochA)>abs(stochB)
                       stochRatio = abs(stochA/stochB);
                       newModel.S(:,rxnParticipation(1)) = newModel.S(:,rxnParticipation(1)) + stochRatio * newModel.S(:,rxnParticipation(2));
                  else
                       stochRatio = abs(stochB/stochA);
                       newModel.S(:,rxnParticipation(1)) = stochRatio * newModel.S(:,rxnParticipation(1)) + newModel.S(:,rxnParticipation(2));
                       fullSolution(rxnParticipation(1)) = fullSolution(rxnParticipation(1))/stochRatio;
                  end
                   newModel.rxns{rxnParticipation(1)} = [newModel.rxns{rxnParticipation(1)} '+' newModel.rxns{rxnParticipation(2)}];

                   newModel = removeIndicatedRxns(newModel, rxnParticipation(2));
                   newModel = removeIndicatedMets(newModel, i);
                   fullSolution(rxnParticipation(2)) = [];
               end
           end
        end
    end

    if rmDuplicates
        [newModel, fullSolution] = removeDuplicates(newModel, fullSolution);
    end
    
    %smallSolution = filterLowFlux(newModel, fullSolution, 0.01);

end



function [model, fullSolution] = removeCompartments(model, fullSolution)
    metaboliteNames = unique(model.metNames);
    for i = 1:length(metaboliteNames)
       matchingMets =  find(ismember(model.metNames, metaboliteNames{i}));
       model.S(matchingMets(1), :) = sum(model.S(matchingMets,:),1);
       model.S(matchingMets(2:end), :) = 0;
    end
    
    emptyMets = sum(abs(model.S),2)==0;
    model = removeIndicatedMets(model, emptyMets);
    
    emptyRxns = sum(abs(model.S),1)==0;
    model = removeIndicatedRxns(model, emptyRxns);
    fullSolution(emptyRxns) = [];  
    
    
end

function [model, fullSolution] = removeDuplicates(model, fullSolution)
    S = model.S;
    storeDuplicates = [];
    for i = length(model.rxns):-1:2
        for j = 1:(i-1)
           if S(:,j) == S(:,i)
               fullSolution(j) = fullSolution(j) + fullSolution(i);             
               storeDuplicates = [storeDuplicates; i];
           end
        end
    end
     model = removeIndicatedRxns(model, storeDuplicates);
    fullSolution(storeDuplicates) = [];   
end

function model = removeIndicatedMets(model, removeMets)   
    %mets
    model.S(removeMets,:) = [];
    model.mets(removeMets) = [];
    model.b(removeMets) = [];
    model.metNames(removeMets) = [];
    model.metComps(removeMets) = [];
    model.metFormulas(removeMets) = [];
    model.metMiriams(removeMets) = [];
    model.inchis(removeMets) = [];
end

function model = removeIndicatedRxns(model, removeRxns)
    %rxns
    model.rxns(removeRxns) = [];
    model.S(:,removeRxns) = [];
    model.lb(removeRxns) = [];
    model.ub(removeRxns) = [];
    model.c(removeRxns) = [];
    model.rev(removeRxns) = [];
    model.rxnNames(removeRxns) = [];
    model.rxnComps(removeRxns) = [];
    model.subSystems(removeRxns) = [];
    model.eccodes(removeRxns) = [];
    model.grRules(removeRxns) = [];
end

function result = areIntegers(A, B)
    result = false;
    if mod(A,1) == 0 && mod(B,1) == 0
        result = true;
    end
end

function [model, fullSolution] = makePositiveModel(model, fullSolution)
    %reverse negative reactions for consistent directionality
    negativeRxns = fullSolution<0;
    model.S(:,negativeRxns) = model.S(:,negativeRxns)*-1;
    fullSolution(negativeRxns) = -fullSolution(negativeRxns);
end