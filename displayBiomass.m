clc
load('model/genericHuman2')
addpath('src')
model = setupBiomass(model, 48, 1);

fprintf('Metabolite\tConcentration\tMol weight\tg/gdw\n');
displayReactionComponents(model, 'HumanGrowth')

