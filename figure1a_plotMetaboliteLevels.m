clc
load('model/genericHuman2')
addpath('src')
global massPerCell
massPerCell = 1000 * 426.8 * 10^-12; %mg per cell

showVitamins = true;

model = setupBiomass(model, 48, 1);

celltype = 'hepg2'; %hepg2 or huh7
condition = '22'; %0, 6 or 22 mM glucose
[expData, volumePoints, growthdat] = loadExpdata('data', celltype, condition);
[fluxMets, fluxValues] = loadFluxes('fluxvalues', celltype, condition);

outputMets = ['biomass[s]';expData.keys'];

halflife = loadHalflife(outputMets);

startT = 0;
tspan = [startT 100];
init = getInitialConditions(expData,growthdat,startT);

haltCrit = zeros(length(outputMets),1);
if strcmp(celltype, 'hepg2')
    haltCrit(findIndex(outputMets, 'glutamine[s]')) = 1200;
else
    haltCrit(findIndex(outputMets, 'glutamine[s]')) = 1400;
end
    
%prevent zero metabolites from triggering growth phase shift
haltCrit(findIndex(outputMets, 'urea[s]')) = -1000;
haltCrit(findIndex(outputMets, '5-oxoproline[s]')) = -1000;
haltCrit(findIndex(outputMets, 'NH3[s]')) = -1000;
haltCrit(findIndex(outputMets, 'asparagine[s]')) = -1000; 
haltCrit(findIndex(outputMets, 'aspartate[s]')) = -1000;
haltCrit(findIndex(outputMets, 'choline[s]')) = -1000;
haltCrit(findIndex(outputMets, 'pantothenate[s]')) = -1000;

newFluxValues = fluxValues;
%newFluxValues = fitMetabolicFluxes(model, tspan, fluxMets, newFluxValues, halflife, outputMets, init, haltCrit, volumePoints, expData);

[t, y, yconc, breakPoints, mu] = fullSimulation(model, tspan, fluxMets, newFluxValues, halflife, outputMets, init, haltCrit, volumePoints, true);

%%
close all
dict = expData.keys;
%dict = {'glucose[s]', 'L-lactate[s]', 'pyruvate[s]', 'glutamine[s]', 'glutamate[s]', 'alanine[s]', 'arginine[s]', 'asparagine[s]', 'aspartate[s]', 'cystine[s]',  'glycine[s]', 'histidine[s]', 'isoleucine[s]', 'leucine[s]', 'lysine[s]', 'methionine[s]', 'phenylalanine[s]', 'serine[s]',  'threonine[s]', 'tyrosine[s]', 'valine[s]'};
dict = {'glucose[s]', 'L-lactate[s]', 'glutamine[s]', 'glutamate[s]', 'alanine[s]', 'phenylalanine[s]'};
%dict = {'glucose[s]', 'L-lactate[s]', 'glutamine[s]', 'alanine[s]'};

numberOfCols = 3;
plotResults(t, y, yconc, breakPoints, mu, expData, growthdat, dict, numberOfCols)

addpath('additionalSamples')
addpath('additionalSamples/src')
%overlayData(dict, condition, numberOfCols)

%getGrowthLimitingMet(model, fluxMets, fluxValues(:,1)/1000)

%%


if showVitamins
    figure()
    dict = {'choline[s]', 'pantothenate[s]', 'folate[s]', 'nicotinamide[s]',  'riboflavin[s]', '5-oxoproline[s]'};
    dict = {'choline[s]', 'pantothenate[s]', 'nicotinamide[s]', '5-oxoproline[s]', 'NH3[s]'};

    numberOfCols = length(dict);
    plotResults(t, y, yconc, breakPoints, mu, expData, growthdat, dict, numberOfCols)
end
