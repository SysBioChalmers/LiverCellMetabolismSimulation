close all
conditions = {
    %'hepg2-0mm-.txt'
    'hepg2-6mm-.txt'
    'hepg2-6mm-.txt' 
    %'huh7-0mm-.txt'
    %'huh7-22mm-.txt'
    };

conditionName = {
    %'hepG2 0mM'
    'hepG2 6mM'    
    %'huh7 0mM'
    'huh7 22mM'
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
    [fluxMets, fluxValues] = loadFluxes('../../fluxvalues', conditions{i});
    fluxMets = strrep(fluxMets, '[s]', '');
    fluxMets{findIndex(fluxMets, 'L-lactate')} = 'lactate';
    plotValues = fluxValues(:,fluxProfile(i));
    for j = 1:length(fluxMets)
        if isKey(result, fluxMets{j})
            tmp = result(fluxMets{j});
        else
            tmp = zeros(1,length(conditions));         
        end
        tmp(i) = plotValues(j);
        result(fluxMets{j}) = tmp;           
    end
end    

allMets = result.keys';
allMets(findIndex(allMets, 'sn-glycerol-3-PC')) = [];
allMets(findIndex(allMets, 'albumin')) = [];

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

[crap, indx] = sort(median(plotData,2), 'descend');


plotData = plotData(indx,:);
allMets = allMets(indx);

largeMets = abs(median(plotData,2))>6*median(abs(plotData(:)));
smallMets = not(largeMets);

    exMap = [67 116 160
             80 137 188
             91 155 213
             151 185 224
             190 209 234]/255;


colormap(exMap)
barh(plotData(largeMets,:), 1, 'LineStyle','none')
legend('Normal', 'Rapid')
%plot(median(plotData(largeMets,:),2),1:sum(largeMets), 'x', 'LineWidth',2,  'color', [0.8500    0.3250    0.0980]);
yticks(1:sum(largeMets))
set(gca,'yticklabel',allMets(largeMets))
ylim([0 sum(largeMets)]+0.5)
xlabel('flux')

xdata = plotData(largeMets,:);

for i = 1:length(xdata)
    ratio = xdata(i,2)./xdata(i,1);
    if sign(xdata(i,2))<0
        side = 50;
    else
        side = -250;
    end
    text(side, i, sprintf('%2.1f', ratio))
end
set(findall(gcf,'-property','FontSize'),'FontSize',15)

figure()
hold all
smallMets = and(smallMets, plotData(:,1)<0);

scatter(plotData(smallMets,1), plotData(smallMets,2),'filled')
%poly = polyfit(plotData(smallMets,1), plotData(smallMets,2), 1);
dlm = fitlm(plotData(smallMets,1),plotData(smallMets,2),'Intercept',false);
slope = predict(dlm,1);

xdata = linspace(-280/slope, 0);
ydata = slope*xdata;


plot(xdata, ydata, '-', 'linewidth', 2, 'color', 0.5 * [1 1 1])
r2 = dlm.Rsquared.Ordinary;
eqn = sprintf('y=%2.2fx\nR2=%2.2f', slope, r2);

text(-40,min(ydata)*0.5, eqn);


metNames = allMets(smallMets);
xpoints = plotData(smallMets,1);
ypoints = plotData(smallMets,2);

for i = 1:length(metNames)
    text(xpoints(i)+2,ypoints(i), metNames{i});
end
xlabel('Normal phase')
ylabel('Rapid phase')
xlim([-50 15])
ylim([-280 0])
set(findall(gcf,'-property','FontSize'),'FontSize',15)
