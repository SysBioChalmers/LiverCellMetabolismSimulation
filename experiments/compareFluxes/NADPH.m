close all
clf
addpath('src')
load('model/genericHuman2')
conditions = {
    'hepg2-0mm-.txt'
    'hepg2-6mm-.txt'
    'hepg2-6mm-.txt'
    %'hepg2-22mm-.txt' 
    %'huh7-0mm-.txt'
    %'huh7-22mm-.txt'
    };

conditionName = {
    'hepG2 0mM'
    'hepG2 6mM'
    'hepG2 6mM (rapid)'
    %'hepG2 22mM'
    %'huh7 0mM'
    %'huh7 22mM'
    };

fluxProfile = [
    3
    2
    1
    %2
    %2
    %1
    ];

flux = zeros(length(model.rxns),length(conditions));

model = setupBiomass(model, 140, 0.5, 0.83);

for i = 1:length(conditions)
    [fluxMets, fluxValues] = loadFluxes('fluxvalues', conditions{i});
    model = bindFBA(model, fluxMets, fluxValues(:,fluxProfile(i))/1000);
    solution = solveLinMin(model,1);
    flux(:,i) = solution.x;
end    


hold on    
subplot(10,1,1:5);
results = zeros(4, length(conditions));
for i = 1:length(conditions)
    [outflux, eqnsOut] = calculateNADPHValues(model, flux(:,i), 'in');
    results(:,i) = outflux;    
end

plotBarWithPercent(results, eqnsOut, conditionName)    
xticks([])
subplot(10,1,6:10);
hold on    

results = zeros(4, length(conditions));
for i = 1:length(conditions)
    [outflux, eqnsOut] = calculateNADPHValues(model, flux(:,i), 'out');
    results(:,i) = outflux;    
end

plotBarWithPercent(results, eqnsOut, conditionName)    

xlabel('NADPH')

