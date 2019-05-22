clc
load('../model/genericHuman2')
addpath('../src')

close all
celltype = {
    'hepG2'
    'hepG2' 
    'hepG2'    
    };

conditions = [
    200
    100
    0
    ];

conditionName = {
    'hepG2 200SSZ'
    'hepG2 100SSZ'
    'hepG2 0SSZ'  
    };

fluxProfile = [
    2
    2
    %2
    2
    %2
    ];

result = containers.Map();
growthRate = zeros(length(conditions),1);

model = setupBiomass(model, 50, 0.5);

for i = 1:length(conditions)
    fileName = ['output/hepg2-' num2str(conditions(i)) 'SSZ.tsv'];
    raw = IO(fileName);
    fluxMets = raw(2:end,1);
    plotValues = cell2nummat(raw(2:end,2));
    plotError = cell2nummat(raw(2:end,3:4));

    fluxMetName = strrep(fluxMets, '[s]', '');
    fluxMetName{findIndex(fluxMetName, 'L-lactate')} = 'lactate';
    
    for j = 1:length(fluxMetName)
        if isKey(result, fluxMetName{j})
            tmp = result(fluxMetName{j});
        else
            tmp = zeros(3,length(conditions));         
        end
        tmp(1,i) = plotValues(j);
        tmp(2,i) = plotError(j,1);
        tmp(3,i) = plotError(j,2);
        result(fluxMetName{j}) = tmp;           
    end
    
    %Set insignificant fluxes to 0
%     aroundZero = sign(fluxError(:,1)) == -sign(fluxError(:,2));
%     plotValues(aroundZero) = 0;
    
    %Calculate growth rate
    model = bindFBA(model, fluxMets, plotValues/1000);
    solution = solveLin(model,1);
    growthRate(i) = -solution.f;
    
end    



growthData = result('biomass');


allMets = result.keys';
allMets(findIndex(allMets, 'sn-glycerol-3-PC')) = [];
allMets(findIndex(allMets, 'albumin')) = [];
allMets(findIndex(allMets, 'biomass')) = [];
%allMets(findIndex(allMets, 'tryptophan')) = [];

plotData = zeros(length(allMets), length(conditions));
plotErrorLB = zeros(length(allMets), length(conditions));
plotErrorUB = zeros(length(allMets), length(conditions));



for i = 1:length(allMets)
    tmp = result(allMets{i});
    plotData(i,:) = tmp(1,:);
    plotErrorLB(i,:) = tmp(2,:);
    plotErrorUB(i,:) = tmp(3,:);
end

%normalize
% factor = abs(median(plotData));
% plotData = plotData./repmat(factor, size(plotData,1),1);


%Reverse data for better plot
plotData = plotData(:,size(plotData,2):-1:1);
plotErrorLB = plotErrorLB(:,size(plotErrorLB,2):-1:1);
plotErrorUB = plotErrorUB(:,size(plotErrorUB,2):-1:1);
growthData = growthData(:,size(growthData,2):-1:1);

conditionName = flipud(conditionName);
growthRate = flipud(growthRate);
[crap, indx] = sort(median(plotData,2), 'descend');


plotData = plotData(indx,:);
plotErrorLB = plotErrorLB(indx,:);
plotErrorUB = plotErrorUB(indx,:);
allMets = allMets(indx);

largeMets = median(abs(plotData),2)>3*median(abs(plotData(:)));
smallMets = not(largeMets);





%%         
hold all
exMap = [%211 111 41
         80 137 188
         %227 120 48
         151 185 224
         190 209 234
         ];
exMap=exMap/255;
secColor = [0.8500    0.3250    0.0980];

xvals = makeErrorBarPlot(growthData(1,:), growthData(2,:), growthData(3,:),'growth', exMap);
plot(xvals, growthRate, '^', 'markerfacecolor', secColor, 'markeredgecolor', secColor)
%bar([growthRate'; zeros(1,length(growthRate))], 1, 'LineStyle','none')
ylabel('specific growth rate 1/h')
set(gca,'FontSize', 15, 'FontWeight', 'bold');
ylim([0 0.08])
xlim([0.5 1.5])
figure

makeErrorBarPlot(plotData(largeMets,:), plotErrorLB(largeMets,:), plotErrorUB(largeMets,:), allMets(largeMets), exMap);
xlim([0 sum(largeMets)]+0.5)
xtickangle(90)
ylabel('Flux')
set(gca,'FontSize', 15, 'FontWeight', 'bold');

figure
hold all
makeErrorBarPlot(plotData(smallMets,:), plotErrorLB(smallMets,:), plotErrorUB(smallMets,:), allMets(smallMets), exMap);

legend(conditionName);
legend boxoff

%plot(median(plotData(smallMets,:),2),1:sum(smallMets), 'x', 'LineWidth',2, 'color', [0.8500    0.3250    0.0980]);
ylabel('Flux')

