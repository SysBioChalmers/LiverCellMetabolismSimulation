clc
clf
load('model/genericHuman2')
addpath('src')
addpath('sensitivity analysis')

[fluxMets, fluxValues] = loadFluxes('fluxvalues', 'hepg2-6mm-.txt');
model = setupBiomass(model, 150, 0.5);
model = bindFBA(model, fluxMets, fluxValues(:,2)/1000);

curSol = solveLinMin(model);
fluxes = curSol.x;



[normalizedGrowth, xValues, reactionNumbers] = ASA(model, fluxes, 10);


%save('sensitivity analysis/sensitivityProfilesMaxFlux.mat','growthRates','reactionNumbers', 'simulationSteps')


%%
hold all
%Displace the lines verticly at random by 1% to prevent overlapping lines
randValues = 0.02 * rand(length(reactionNumbers), 1)-0.01; 
normalizedGrowthAndDisplacement = normalizedGrowth + repmat(randValues, 1, length(xValues));

for i = 1:size(normalizedGrowthAndDisplacement,1)
    plot(xValues, normalizedGrowthAndDisplacement(i,:), 'color', [0.3 0.3 0.3 0.1]);
end

colors = [
         0    0.4470    0.7410
    0.8500    0.3250    0.0980
    0.9290    0.6940    0.1250
    0.4940    0.1840    0.5560
    0.4660    0.6740    0.1880
    0.3010    0.7450    0.9330
    0.6350    0.0780    0.1840];

conditions = {'glucose[s]', 'glutamine[s]', 'O2[s]'};
interestingReactions = getBounds(model, conditions);

%plot([0 1], [0, 1], '--', 'color', [0 0 0], 'linewidth', 2);


legendIds = zeros(length(interestingReactions),1);
for i = 1:length(interestingReactions)
    curRxn = ismember(reactionNumbers, interestingReactions(i));
    legendIds(i) = plot(xValues, normalizedGrowth(curRxn,:), 'linewidth', 3, 'color', colors(i,:));
end
legend(legendIds,conditions);
legend boxoff


axis equal
xlim([0 1])

xlabel('restriction')
ylabel('growth')
ylim([0 inf])

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
conditions = {'glutamate[s]'};
interestingReactions = getBounds(model, conditions);
curRxn = ismember(reactionNumbers, interestingReactions);
plot(xValues, normalizedGrowth(curRxn,:), 'linewidth', 3, 'color', colors(i,:))

plot([0 1], [0,1], 'k--')
legend boxoff

%%
rxnIds = model.rxns(reactionNumbers);
xLabels = cellstr(num2str(round(xValues',2)));

interestingRxns = 1:length(normalizedGrowth);
nonLinearDecline = normalizedGrowth(:,3)>1.001*xValues(:,3);
nonZeroDecline = normalizedGrowth(:,1)<0.9;
interestingRxns = and(nonLinearDecline,nonZeroDecline);

h = clustergram(normalizedGrowth(interestingRxns,:), 'Cluster', 'column', 'ColumnLabels', xLabels, 'RowLabels',rxnIds(interestingRxns), 'Colormap', 'redbluecmap')


