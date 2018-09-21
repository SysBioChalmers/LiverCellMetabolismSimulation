function [newModel, fullSolution] = simplifyFluxResult(model, fullSolution, rmComps)
% simplifyFluxResult
% Simplifies the model by joining linear pathways to lumped reactions
%
%   model               a model struct
%   fullSolution        a vector with fluxes
%   rmComps             boolean, remove compartments
%   
%   Avlant Nilsson, 2017-05-17
%
fluxThresh = 10^-6;

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
                   
               newModel = removeIndicatedRxns(newModel, rxnParticipation(2));
               newModel = removeIndicatedMets(newModel, i);
               fullSolution(rxnParticipation(2)) = [];
           end
       end
    end
end 

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

% function [model, fullSolution] = removeDuplicates(model, fullSolution)
%     absS = abs(model.S);
%     storeDuplicates = [];
%     for i = length(model.rxns):-1:2
%         for j = 1:(i-1)
%            if absS(:,j) == absS(:,i)
%                metabolites = find(absS(:,i));
%                absVal = model.S(metabolites,j)./model.S(metabolites,i);
%                direction = sign(sum(absVal));
%                fullSolution(j) = fullSolution(j) + direction*fullSolution(i);
%                
%                storeDuplicates = [storeDuplicates; i];
%            end
%         end
%     end
%      model = removeIndicatedRxns(model, storeDuplicates);
%     fullSolution(storeDuplicates) = [];   
% end

function model = removeIndicatedMets(model, removeMets)   
    %mets
    model.S(removeMets,:) = [];
    model.mets(removeMets) = [];
    model.b(removeMets) = [];
    model.metNames(removeMets) = [];
    model.metComps(removeMets) = [];
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