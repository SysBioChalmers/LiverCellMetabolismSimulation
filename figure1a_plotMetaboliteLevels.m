clc
load('model/genericHuman2')
addpath('src')
global massPerCell
massPerCell = 1000 * 426.8 * 10^-12; %mg per cell


model = setupBiomass(model, 150, 0.5);



celltype = 'hepg2';
condition = '6mm';
[expData, volumePoints, growthdat] = loadExpdata('data', celltype, condition);
[fluxMets, fluxValues] = loadFluxes('fluxvalues', [celltype '-' condition '-.txt']);

outputMets = ['biomass[s]';expData.keys'];

halflife = loadHalflife(outputMets);

startT = 0;
tspan = [startT 100];
init = getInitialConditions(expData,growthdat,startT);

haltCrit = zeros(length(outputMets),1);


haltCrit(findIndex(outputMets, 'glutamine[s]')) = 1200;

%prevent zero metabolites from triggering
haltCrit(findIndex(outputMets, 'urea[s]')) = -1000;
haltCrit(findIndex(outputMets, '5-oxoproline[s]')) = -1000;
haltCrit(findIndex(outputMets, 'NH3[s]')) = -1000;
haltCrit(findIndex(outputMets, 'asparagine[s]')) = -1000; 
haltCrit(findIndex(outputMets, 'aspartate[s]')) = -1000;
%haltCrit(findIndex(outputMets, 'choline[s]')) = -1000;


newFluxValues = fluxValues;
%newFluxValues = fitMetabolicFluxes(model, tspan, fluxMets, newFluxValues, halflife, outputMets, init, haltCrit, volumePoints, expData);

[t, y, yconc, breakPoints, mu] = fullSimulation(model, tspan, fluxMets, newFluxValues, halflife, outputMets, init, haltCrit, volumePoints, true);


close all
dict = expData.keys;
dict = {'glucose[s]', 'L-lactate[s]', 'pyruvate[s]', 'glutamine[s]', 'glutamate[s]', 'alanine[s]', 'arginine[s]', 'asparagine[s]', 'aspartate[s]', 'cystine[s]',  'glycine[s]', 'histidine[s]', 'isoleucine[s]', 'leucine[s]', 'lysine[s]', 'methionine[s]', 'phenylalanine[s]', 'serine[s]',  'threonine[s]', 'tyrosine[s]', 'valine[s]'};
%dict = {'glucose[s]', 'L-lactate[s]', 'pyruvate[s]', 'glutamine[s]', 'glutamate[s]', 'alanine[s]'};
%dict = { 'arginine[s]', 'asparagine[s]', 'aspartate[s]', 'taurine[s]', 'cystine[s]',  'glycine[s]', 'histidine[s]', 'isoleucine[s]', 'leucine[s]', 'lysine[s]', 'methionine[s]', 'phenylalanine[s]', 'serine[s]',  'threonine[s]', 'tyrosine[s]', 'valine[s]'};
%dict = {'choline[s]', 'pantothenate[s]', 'folate[s]', 'nicotinamide[s]',  'riboflavin[s]'};
%dict = {'taurine[s]'};
dict = {'glucose[s]', 'L-lactate[s]', 'glutamine[s]', 'pyruvate[s]', 'glutamate[s]', 'alanine[s]'};

plotResults(t, y, yconc, breakPoints, mu, expData, growthdat, dict)

%getGrowthLimitingMet(model, fluxMets, fluxValues(:,1)/1000)