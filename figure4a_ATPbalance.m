close all
clf
addpath('src')
load('model/genericHuman2')
conditions = {
    'hepg2-6mm-.txt'
    'hepg2-6mm-.txt'
    %'hepg2-22mm-.txt' 
    %'huh7-22mm-.txt'
    'hepg2-0mm-.txt'    
    %'huh7-0mm-.txt'
    };

conditionName = {
    'hepG2 6mM (rapid)'
    'hepG2 6mM'
    %'hepG2 22mM'
    %'huh7 6mM & 22mM'
    'hepG2 0mM'  
    %'huh7 0mM'
    };

fluxProfile = [
    1
    2
    %2
    %2
    2
    %2
    ];

model = setupBiomass(model, 150, 0.5);


flux = zeros(length(model.rxns),length(conditions));


for i = 1:length(conditions)
    [fluxMets, fluxValues] = loadFluxes('fluxvalues', conditions{i});
    model = bindFBA(model, fluxMets, fluxValues(:,fluxProfile(i))/1000);
    solution = solveLinMin(model,1);
    flux(:,i) = solution.x;
end    


%%
subplot(10,1,1:5);
hold on 
set(gca,'DefaultTextFontSize',12)
results = zeros(4, length(conditions));
for i = 1:length(conditions)
    [outflux, eqnsOut] = calculateATPValues(model, flux(:,i), 'in');
    results(:,i) = outflux;    
end

plotBarWithPercent(results, eqnsOut, conditionName)    
xticks([])
set(findall(gcf,'-property','FontSize'),'FontSize',15)

pos = get(gca, 'Position');
pos(1) = pos(1)*2.5;
pos(3) = pos(3)*0.8;
set(gca, 'Position', pos)

subplot(10,1,6:10);
hold on    


results = zeros(5, length(conditions));
for i = 1:length(conditions)
    [outflux, eqnsOut] = calculateATPValues(model, flux(:,i), 'out');
    results(:,i) = outflux;    
end

plotBarWithPercent(results, eqnsOut, conditionName)    

pos = get(gca, 'Position');
pos(1) = pos(1)*2.5;
pos(3) = pos(3)*0.8;
set(gca, 'Position', pos)

estimatedProteinTurnover = 2.57;
plot(estimatedProteinTurnover * [1 1], [0 4], 'k-')

xlabel('ATP')
set(findall(gcf,'-property','FontSize'),'FontSize',15)
