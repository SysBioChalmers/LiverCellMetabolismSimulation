clc
clf
load('model/genericHuman2')
addpath('src')
addpath('sensitivity analysis')


resolution = 40;
celltype = 'hepg2';
fluxProfile = 2;
condition = 22;
tresh = 10^-6;

if and(strcmp(celltype, 'hepg2'), fluxProfile == 2)
    fileName = ['confidenceIntervalls\output\hepg2-' num2str(condition) '.tsv'];
    raw = IO(fileName);
    fluxMets = raw(2:end,1);
    fluxIn = cell2nummat(raw(2:end,2))/1000;
else
    [fluxMets, fluxValues] = loadFluxes('fluxvalues', celltype, condition);
    fluxIn = fluxValues(:,fluxProfile)/1000;
end

model = setupBiomass(model, 48, 1);
model = bindFBA(model, fluxMets, fluxIn);

curSol = solveLinMin(model);
fluxes = curSol.x;
growthRate = -curSol.f;

reactionNumbers = abs(fluxes)>tresh;
transpReactions = and(reactionNumbers, contains(model.subSystems, 'Transport'));

[normalizedGrowth, xValues, reactionNumbers] = ASA(model, fluxes, resolution, true, reactionNumbers, find(transpReactions));

hyperSensitive = normalizedGrowth(:,2) < xValues(2)*0.9;
%constructEquations(model, reactionNumbers(hyperSensitive))

%save('sensitivity analysis/sensitivityProfilesMaxFlux.mat','growthRates','reactionNumbers', 'simulationSteps')


%%
hold all
set(gca,'DefaultTextFontSize',14)
set(gca,'DefaultTextFontWeight','bold')

%Displace the lines verticly at random by 0.5% to prevent overlapping lines
randValues = 0.01 * rand(length(reactionNumbers), 1)-0.005; 
normalizedGrowthAndDisplacement = normalizedGrowth + repmat(randValues, 1, length(xValues));

for i = 1:size(normalizedGrowthAndDisplacement,1)
    plot(flip(xValues), growthRate*normalizedGrowthAndDisplacement(i,:), 'color', [0.3 0.3 0.3 0.1]);
end

colors = [
         0    0.4470    0.7410
    0.8500    0.3250    0.0980
    0.9290    0.6940    0.1250
    0.4940    0.1840    0.5560
    0.4660    0.6740    0.1880
    0.3010    0.7450    0.9330
    0.6350    0.0780    0.1840];

condition = {'ATP synthase', 'Hexokinase', 'Lactate dehydrogenase'};
producersOfATP = {
    'HMR_6916'
    'HMR_4394'
    'HMR_4388'
    };

legendIds = zeros(length(producersOfATP),1);
for i = 1:length(producersOfATP)
    curNr = findIndex(model.rxns, producersOfATP{i});
    curRxn = findIndex(reactionNumbers, curNr);
    legendIds(i) = plot(flip(xValues), growthRate*normalizedGrowth(curRxn,:), 'linewidth', 3, 'color', colors(i,:));
end
legend(legendIds,condition, 'location', 'SW');
legend boxoff


%experimentalGrowh = [0.0251 0.0156]; %No glucose, %No glutamine 

%experimental glucose
%plot(1, experimentalGrowh(1), '*', 'linewidth', 3, 'color', colors(1,:));

%experimental glutamine
%plot(1, experimentalGrowh(2), '*', 'linewidth', 3, 'color', colors(2,:));



%axis equal
xlim([0 1])

xlabel('reaction inhibition')
ylabel('specific growth rate [h-1]')
ylim([0 growthRate])

%%
figure()
hold all
BCAA = {'HMR_3777' %isoleucine
        'HMR_6923' %leucine
        'HMR_3747'}; %valine
    
labels = {'BCAT1 (isoleucine)', 'BCAT1 (leucine)', 'BCAT1 (valine)'};    
    
for i = 1:length(BCAA)
    curRxn = findIndex(model.rxns, BCAA{i});
    curRxn = ismember(reactionNumbers, curRxn);
    plot(xValues, normalizedGrowth(curRxn,:), 'linewidth', 3, 'color', colors(i,:))
end
plot([0 1], [0,1], 'k--')
legend(labels, 'location', 'se')
legend boxoff


%%
figure()
hold all
condition = {'glutamate[s]'};
interestingReactions = getBounds(model, condition);
curRxn = ismember(reactionNumbers, interestingReactions);
plot(xValues, normalizedGrowth(curRxn,:), 'linewidth', 3)

plot([0 1], [0,1], 'k--')
legend boxoff

%%
rxnIds = model.rxns(reactionNumbers);
%rxnIds = constructEquations(model, reactionNumbers);
xLabels = cellstr(num2str(round(xValues',2)));

interestingRxns = 1:length(normalizedGrowth);
% nonLinearDecline = normalizedGrowth(:,3)>1.001*xValues(:,3);
% nonZeroDecline = normalizedGrowth(:,1)<0.9;
% interestingRxns = and(nonLinearDecline,nonZeroDecline);

h = clustergram(normalizedGrowth(interestingRxns,:)', 'Cluster', 'row', 'RowLabels', xLabels, 'ColumnLabels',rxnIds(interestingRxns), 'Colormap', 'redbluecmap')


