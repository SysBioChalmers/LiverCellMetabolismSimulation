clc
load('model/genericHuman2')
addpath('src')
model = setupBiomass(model, 150, 0.5);

fprintf('Metabolite\tConcentration\tMol weight\tg/gdw\n');
displayReactionComponents(model, 'HumanGrowth')

