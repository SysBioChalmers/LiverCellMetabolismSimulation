close all
raw = IO('SSZ/growthSSZ.txt');
plotColors = get(gca,'ColorOrder');

columns = raw(1,:);
data = cell2nummat(raw(2:end,:));
data(isnan(data(:,3)),:) = []; %remove missing values

%Remove duplicated timepoints at baseline
data(and(data(:,1) == 0, data(:,2)>0),:) = [];

%Remove points after glucose depletion
data(data(:,1)>50,:) = [];

time = data(:,1);
condition = data(:,2);
cells = data(:,3);

hold all
logCells = log(cells);
dose = condition/1000; %unit mol
data = table(time, logCells, dose);
mdlL = fitlm(data, 'logCells~time:(1+dose)');
mdlL

intercept = mdlL.Coefficients.Estimate(1);
slope = mdlL.Coefficients.Estimate(2);
doseEffect = mdlL.Coefficients.Estimate(3);
SE = mdlL.Coefficients.SE;
pvalue = mdlL.Coefficients.pValue(3);

allConditons = unique(condition);

tpoints = linspace(0, 50);

for i = 1:length(allConditons)
    filter = condition == allConditons(i);
    dataX = time(filter);
    dataY = logCells(filter);
    scatter(dataX, dataY, 'filled', 'markerfacecolor', plotColors(i,:), 'HandleVisibility','off')
    
    curMu = slope + doseEffect * allConditons(i)/1000;
    ypoints = intercept + curMu*tpoints;
    %ypoints = exp(intercept) * exp(curMu*tpoints);
    plot(tpoints, ypoints, 'color', plotColors(i,:), 'linewidth', 2);
end

legend({'SSZ 0 mM', 'SSZ 0.1 mM', 'SSZ 0.2 mM'}, 'location', 'nw')
legend boxoff
xlabel('time [h]')
ylabel('log(cells)')
xlim([0 50])

text(20, 13, sprintf('k = %2.3f + (%2.3f x dose)\ny = %2.2f + kt\np = %2.2e', slope, doseEffect, intercept, pvalue))



