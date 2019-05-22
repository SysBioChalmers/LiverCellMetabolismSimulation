metData = IO('SSZ/SSZdata.tsv');
colNames = metData(1,:);

includeMets = {'Cystine' 'Glutamic acid' 'Alanine' 'Arginine' 'Aspartic acid' 'Glutamine' 'Glycine' 'Histidine' 'Iso-leucine' 'Leucine' 'Lysine' 'Methionine' 'Ornithine' 'Phenylalanine' 'Proline' 'Serine' 'Threonine' 'Tryptophan' 'Tyrosine' 'Valine'};

primColor = [80 137 188]/255;
secColor = [0.8500    0.3250    0.0980];


% exMap = [%211 111 41
%          80 137 188
%          %227 120 48
%          151 185 224
%          190 209 234
%          ];

%Remove mets
metNames = metData(:,ismember(colNames,'met'));
metData(not(ismember(metNames, includeMets)),:) = [];

%Extract condition
con = cell2nummat(metData(:,ismember(colNames,'condition')))/1000; %mM

%Extract timepoints
time = cell2nummat(metData(:,ismember(colNames,'time')));

%Extract values
value = cell2nummat(metData(:,ismember(colNames,'value')));

%Make met ids
metNames = metData(:,ismember(colNames,'met'));
metId = zeros(length(metNames),1);
for i = 1:length(metId)
    metId(i) = findIndex(includeMets, metNames{i});
end

dataY = value;
dataX = [time con metId];

filter = ismember(dataX(:,1), [0 24.75 45.25]);

dataX = dataX(filter,:);
dataY = dataY(filter);

result = zeros(length(includeMets),4);
meanValues = zeros(length(includeMets),3);
for i = 1:length(includeMets)
    filter = dataX(:,3) == i; 
    time = dataX(filter,1);
    value = dataY(filter);
    dose = dataX(filter,2);
    data = table(time, value, dose);
    mdlL = fitlm(data, 'value~time:(1+dose)');
    %mdlL

    result(i,1) = mdlL.Coefficients.Estimate(1);
    result(i,2) = mdlL.Coefficients.Estimate(2);
    result(i,3) = mdlL.Coefficients.Estimate(3);
    result(i,4) = mdlL.Coefficients.pValue(3);
    
    meanValues(i, 1) = mean(value(and(time==45.25, dose == 0)));
    meanValues(i, 2) = mean(value(and(time==45.25, dose == 100)));
    meanValues(i, 3) = mean(value(and(time==45.25, dose == 200)));
end
pval = result(:,4);

%%

% 
% subplot(1,2,1)
% hold all
% curMet = findIndex(includeMets, 'Cystine');
% filter = and(dataX(:,3) == curMet, dataX(:,1) == 45.25);
% value = dataY(filter);
% dose = dataX(filter,2);
% bar([0 100 200], meanValues(1,:), 'LineStyle','none', 'FaceColor', primColor)
% translucentScatter(dose, value, secColor, 0.8, 1)
% ylim([0 350])
% xlim([-75 275])
% 
% subplot(1,2,2)
% hold all
% curMet = findIndex(includeMets, 'Glutamic acid');
% filter = and(dataX(:,3) == curMet, dataX(:,1) == 45.25);
% value = dataY(filter);
% dose = dataX(filter,2);
% bar([0 100 200], meanValues(2,:), 'LineStyle','none', 'FaceColor', primColor)
% translucentScatter(dose, value, secColor, 0.8, 1)
% ylim([0 350])
% xlim([-75 275])

%%
figure()
plotColors = get(gca,'ColorOrder');
conditions = [0 0.1 0.2];

tpoints = linspace(0, 50);

subplot(1,2,1)
hold all
curMet = findIndex(includeMets, 'Cystine');
filter = dataX(:,3) == curMet;
time = dataX(filter,1);
value = dataY(filter);
dose = dataX(filter,2);
coef = result(curMet,:);

for i = 1:length(conditions)
    filter = (dose == conditions(i));
    scatter(time(filter), value(filter), 'fill', 'markerfacecolor', plotColors(i,:), 'HandleVisibility','off')
    ypoints = coef(1) + (coef(2) + coef(3) * conditions(i)) * tpoints;
    plot(tpoints, ypoints, 'color', plotColors(i,:), 'linewidth', 2);    
end
ylim([0 200])
xlim([0 50])
xlabel('time [h]')
ylabel('concentration [µmol/l]')
text(2, 50, sprintf('k = %2.2f + (%2.2f x dose)\ny = %2.2f + kt\np = %2.2e', coef(2), coef(3), coef(1), coef(4)))
title('Cystine')

subplot(1,2,2)
hold all
curMet = findIndex(includeMets, 'Glutamic acid');
filter = dataX(:,3) == curMet;
time = dataX(filter,1);
value = dataY(filter);
dose = dataX(filter,2);
coef = result(curMet,:);

for i = 1:length(conditions)
    filter = (dose == conditions(i));
    scatter(time(filter), value(filter), 'fill', 'markerfacecolor', plotColors(i,:), 'HandleVisibility','off');    
    ypoints = coef(1) + (coef(2) + coef(3) * conditions(i)) * tpoints;
    plot(tpoints, ypoints, 'color', plotColors(i,:), 'linewidth', 2);
end
ylim([0 400])
xlim([0 50])
xlabel('time [h]')
ylabel('concentration [µmol/l]')
text(2, 300, sprintf('k = %2.2f + (%2.2f x dose)\ny = %2.2f + kt\np = %2.2e', coef(2), coef(3), coef(1), coef(4)))
title('Glutamic acid')

%%
figure()

for i = 1:length(includeMets)
    subplot(4,6,i)
    hold all
    filter = dataX(:,3) == i;
    time = dataX(filter,1);
    value = dataY(filter);
    dose = dataX(filter,2);
    coef = result(i,:);

    for j = 1:length(conditions)
        filter = (dose == conditions(j));
        scatter(time(filter), value(filter), 'fill', 'markerfacecolor', plotColors(j,:), 'HandleVisibility','off');    
        ypoints = coef(1) + (coef(2) + coef(3) * conditions(j)) * tpoints;
        plot(tpoints, ypoints, 'color', plotColors(j,:), 'linewidth', 2);
    end
    xlim([0 50])
    tmp = ylim;
    ylim([0 tmp(2)])
    xlabel('time [h]')
    ylabel('concentration [µmol/l]')
    text(2, 2, sprintf('k = %2.2f + (%2.2f x dose)\ny = %2.2f + kt\np = %2.2e', coef(2), coef(3), coef(1), coef(4)))
    title(includeMets{i})
end