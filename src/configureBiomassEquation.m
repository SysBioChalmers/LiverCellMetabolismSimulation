function model = configureBiomassEquation(model, biomassReaction, fatRatio, growthRelMaintenance)
% configureBiomassEquation
% Sets upp the biomass equation with the calculated growth related
% maintenance and the fat ratio
%
%   model                a model structure
%   biomassReaction      a string reference to the biomass rxn, e.g. 
%                        'human_biomass'
%   fatRatio             a value between 0 and 1 with the ratio
%                        fatmass/total mass
%   growthRelMaintenance The growth related energy expenditure in mMol
%                        ATP/gram dry weight
%   model                an updated model structure
%
%   Avlant Nilsson, 2016-05-16
%
    
    %Identify and clear biomass rxn
    rxnNR = findIndex(model.rxns, biomassReaction);
    model.S(:, rxnNR) = 0;
    
    %Set the individual components of the biomass equation
    model = configureSMatrix(model, fatRatio, biomassReaction, 'biomassFat[c]');    
    model = configureSMatrix(model, 1 - fatRatio, biomassReaction, 'biomassLean[c]');
    model = configureSMatrix(model,  growthRelMaintenance, biomassReaction, 'human_growthMaintainance[c]');
    model = configureSMatrix(model, -1, biomassReaction, 'human_biomass[c]');
    
    
    
end

