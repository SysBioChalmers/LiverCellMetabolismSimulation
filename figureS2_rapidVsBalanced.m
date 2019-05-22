close all
celltype = {
    %'hepg2'
    'hepg2'
    'hepg2' 
    %'huh7'
    %'huh7'
    };

conditions = [
    %'0'
    22
    22 
    %'0'
    %'22'
    ];

conditionName = {
    %'hepG2 0mM'
    'hepG2 22mM (rapid)'    
    'hepG2 22mM (balanced)'    
    %'huh7 0mM'
    };

fluxProfile = [
    %2
    1
    2
    %2
    %1
    ];

result = containers.Map();

for i = 1:length(conditions)
    if and(strcmp(celltype{i}, 'hepg2'), fluxProfile(i) == 2)
        fileName = ['confidenceIntervalls\output\hepg2-' num2str(conditions(i)) '.tsv'];
        raw = IO(fileName);
        fluxMets = raw(2:end,1);
        plotValues = cell2nummat(raw(2:end,2));
        plotError = cell2nummat(raw(2:end,3:4));
    else
        [fluxMets, fluxValues] = loadFluxes('fluxvalues', celltype{i}, conditions(i));
        plotValues = fluxValues(:,fluxProfile(i));
        plotError = nan(length(plotValues),2);
    end
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
end    


allMets = result.keys';
allMets(findIndex(allMets, 'albumin')) = [];

%The absolut value is very small, could be noise
allMets(findIndex(allMets, 'aspartate')) = [];



plotData = zeros(length(allMets), 2);
plotErrorLB = zeros(length(allMets), 2);
plotErrorUB = zeros(length(allMets), 2);

for i = 1:length(allMets)
    tmp = result(allMets{i});
    plotData(i,:) = tmp(1,:);
    plotErrorLB(i,:) = tmp(2,:);
    plotErrorUB(i,:) = tmp(3,:);
end

[crap, indx] = sort(median(plotData,2), 'descend');

plotData = plotData(indx,:);
plotErrorLB = plotErrorLB(indx,:);
plotErrorUB = plotErrorUB(indx,:);
allMets = allMets(indx);


plotData(plotData==0) = 10^-6;

largeMets = log(plotData(:,1)./plotData(:,2))>1;

%%
exMap = [%211 111 41
         80 137 188
         %227 120 48
         151 185 224
         190 209 234
         ];
exMap=exMap/255;
secColor = [0.8500    0.3250    0.0980];
xval = makeErrorBarPlot(plotData(largeMets,:), plotErrorLB(largeMets,:), plotErrorUB(largeMets,:), allMets(largeMets), exMap);

ratios = plotData(largeMets,1)./plotData(largeMets,2);
signOfFlux = sign(plotData(largeMets,1));

for i = 1:length(xval)
    if signOfFlux(i) == 1
       y = -50;
    else
       y = 50;
    end    
    text(xval(i,1), y, sprintf('%2.1f',ratios(i)));
end
ylim([-500 800])
legend(conditionName, 'location', 'ne')
legend boxoff
