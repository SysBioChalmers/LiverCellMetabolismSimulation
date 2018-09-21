clc
load('model/genericHuman2')
addpath('src')

close all
conditions = {
    %'hepg2-6mm-.txt'
    %'hepg2-6mm-.txt'
    'hepg2-6mm-.txt'
    'hepg2-22mm-.txt' 
    'huh7-22mm-.txt'
    'hepg2-0mm-.txt'    
    'huh7-0mm-.txt'
    };

conditionName = {
    'hepG2 6mM'
    'hepG2 22mM'
    'huh7 22mM'
    'hepG2 0mM'  
    'huh7 0mM'
    };

fluxProfile = [
    2
    2
    2
    2
    2
    ];

result = containers.Map();
growthRate = zeros(length(conditions),1);
for i = 1:length(conditions)
    [fluxMets, fluxValues] = loadFluxes('fluxvalues', conditions{i});
    fluxMetName = strrep(fluxMets, '[s]', '');
    fluxMetName{findIndex(fluxMetName, 'L-lactate')} = 'lactate';
    plotValues = fluxValues(:,fluxProfile(i));
    for j = 1:length(fluxMetName)
        if isKey(result, fluxMetName{j})
            tmp = result(fluxMetName{j});
        else
            tmp = zeros(1,length(conditions));         
        end
        tmp(i) = plotValues(j);
        result(fluxMetName{j}) = tmp;           
    end
    
    %Calculate growth rate
    model = bindFBA(model, fluxMets, fluxValues(:,fluxProfile(i))/1000);
    solution = solveLin(model,1);
    growthRate(i) = -solution.f;
    
end    

model = setupBiomass(model, 150, 0.5);

allMets = result.keys';
allMets(findIndex(allMets, 'sn-glycerol-3-PC')) = [];
allMets(findIndex(allMets, 'albumin')) = [];
allMets(findIndex(allMets, 'tryptophan')) = [];

plotData = zeros(length(allMets), length(conditions));
for i = 1:length(allMets)
    plotData(i,:) = result(allMets{i});
end

%normalize
% factor = abs(median(plotData));
% plotData = plotData./repmat(factor, size(plotData,1),1);


%Reverse data for better plot
plotData = plotData(:,size(plotData,2):-1:1);
conditionName = flipud(conditionName);
growthRate = flipud(growthRate);
[crap, indx] = sort(median(plotData,2), 'descend');


plotData = plotData(indx,:);
allMets = allMets(indx);

largeMets = median(abs(plotData),2)>3*median(abs(plotData(:)));
smallMets = not(largeMets);

exMap = [211 111 41
         80 137 188
         227 120 48
         151 185 224
         190 209 234]/255;

        





%%         
hold all
colormap(exMap)
bar([growthRate'; zeros(1,length(growthRate))], 1, 'LineStyle','none')
xlim([0 1]+0.5)
ylabel('specific growth rate 1/h')
set(gca,'FontSize', 15, 'FontWeight', 'bold');
figure

colormap(exMap)
bar(plotData(largeMets,:), 1, 'LineStyle','none')
%plot(median(plotData(largeMets,:),2),1:sum(largeMets), 'x', 'LineWidth',2,  'color', [0.8500    0.3250    0.0980]);
xticks(1:sum(largeMets))
set(gca,'xticklabel',allMets(largeMets))
xlim([0 sum(largeMets)]+0.5)
xtickangle(90)
ylabel('Flux')
set(gca,'FontSize', 15, 'FontWeight', 'bold');
    
figure
hold all
colormap(exMap)
h = bar(plotData(smallMets,:), 1, 'LineStyle','none');
legend(conditionName);

%plot(median(plotData(smallMets,:),2),1:sum(smallMets), 'x', 'LineWidth',2, 'color', [0.8500    0.3250    0.0980]);
legend boxoff
xticks(1:sum(smallMets))
set(gca,'xticklabel', allMets(smallMets))
ylabel('Flux')
xtickangle(90)
set(gca,'FontSize', 15, 'FontWeight', 'bold');
xlim([0 sum(smallMets)]+0.5)
ylim([-35.1 50.1])