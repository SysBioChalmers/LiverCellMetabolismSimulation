clc
close all
addpath('../confidenceIntervalls/src')
addpath('src')
condition = '0';

confidencebound = 0.95; %2 std

massPerCell = 426.8 * 10^-12; %g per cell
massPerMCell = massPerCell*10^6;

color = [93 155 211
         215 86 40
         238 178 32]/256;

%leaving out  'Aspartic acid' 'Asparagine' 'Butyrate', 'Taurine', 'Citrulline' 
metabolitesToPlot = {'Glucose' 'Lactate' 'Pyruvate' 'Alanine' 'Arginine'  'Cystine' 'Glutamic acid' 'Glutamine' 'Glycine' 'Histidine' 'Iso-leucine' 'Leucine' 'Lysine' 'Methionine' 'Ornithine' 'Phenylalanine' 'Proline' 'Serine' 'Threonine' 'Tryptophan' 'Tyrosine' 'Valine'};
%metabolitesToPlot = {'Alanine'};

%No sample removal in this experiment
volData = [24.75    2
           45.25    2
           68.75    2];

startTime = volData(1,1);
volData(:,1) = volData(:,1) - startTime;

addConstraints = NaN(length(metabolitesToPlot)+1,3);
[dataX, dataY, metOrder] = makeDataStructureSSZ(condition, metabolitesToPlot);

%Estimating and constraining STD of cells using only growth data
disp('Estimating STD of growth...')
[dataX, dataY] = makeDataStructureSSZ(condition, {});
tvals = unique(dataX(:,1))';
x0 = estimateInitialX(dataX, dataY, tvals);
modelfun = @(x,y) fitFunction(x, y, tvals, volData(:,2));
[beta, fval] = mLE(dataX, dataY, modelfun, x0, NaN(1,3));
addConstraints(1,3) = beta(1,3);

%Estimating and constraining STD for each metabolite in turn using the data
disp('Estimating STD of each metabolite...')
for i = 2:length(metOrder)
    [dataX, dataY] = makeDataStructureSSZ(condition, metOrder(i));
    tvals = unique(dataX(:,1))';
    x0 = estimateInitialX(dataX, dataY, tvals);
    modelfun = @(x,y) fitFunction(x, y, tvals, volData(:,2));
    curConstraints = NaN(2,3);
    curConstraints(1,:) = addConstraints(1,:);
    [beta, fval] = mLE(dataX, dataY, modelfun, x0, curConstraints);
    addConstraints(i,3) = beta(2,3);
end

%%
allBeta = zeros(length(metabolitesToPlot),3);
allGrowthBeta = zeros(length(metabolitesToPlot),3);

allConf = zeros(length(metabolitesToPlot),2);
allGrowthConf = zeros(length(metabolitesToPlot),2);

disp('Fitting individual metabolites and calculating CIs...')
for i = 1:length(metabolitesToPlot)
    disp(metabolitesToPlot{i})
    curMet = metabolitesToPlot{i};
    curMetNr = findIndex(metOrder, curMet);
    [dataX, dataY, mets] = makeDataStructureSSZ(condition, {curMet});
    tvals = unique(dataX(:,1))';
    x0 = estimateInitialX(dataX, dataY, tvals);
    
    curConstraints = addConstraints([1; curMetNr],:);
    
    modelfun = @(x,y) fitFunction(x, y, tvals, volData(:,2));
    [curBeta, conf] = mLEConfidence(dataX, dataY, modelfun, x0, confidencebound, curConstraints);
    
    allGrowthBeta(i,:) = curBeta(1,:);
    allBeta(i,:) = curBeta(2,:);
    allGrowthConf(i,:) = conf(1,:);
    allConf(i,:) = conf(2,:);
end

%%
figure()
trange = linspace(tvals(1),tvals(end));

for i = 1:length(metabolitesToPlot)
    subplot(5,6,i)
    hold all
    curBeta = [allGrowthBeta(i,:);allBeta(i,:)];
    metName = metabolitesToPlot{i};
    
    [dataX, dataY, metlabels] = makeDataStructureSSZ(condition, metName);
    curMet = findIndex(metlabels, metName);
    filter = dataX(:,2)==curMet;
    xrange = [trange' curMet * ones(length(trange),1)];
    
    ypred = modelfun(curBeta, xrange);
    delta = curBeta(curMet,3) * 1.96;
    
    lower = ypred - delta;
    upper = ypred + delta;
    ul = [lower upper];

    fillArea(startTime + trange, ul, color(1,:))
    scatter(startTime + dataX(filter,1), dataY(filter), 'MarkerFaceColor', color(2,:), 'MarkerEdgeColor', color(2,:),'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 1)
    plot(startTime+ trange, ypred, 'k-', 'linewidth', 2)

    title(metName)
    maxval = ylim;
    ylim([0 maxval(2)])    
end

%%
figure()
hold all

predFlux = allGrowthBeta(:,2);
confFlux = allGrowthConf;

xLocation = 1:length(predFlux);
bar(xLocation,predFlux)
errorbar(xLocation, predFlux, confFlux(:,1)-predFlux, confFlux(:,2)-predFlux,'k.');

jointGrowth = mean(predFlux);
jointConfidence = mean(confFlux,1);

bar(length(xLocation)+1, jointGrowth, 'facecolor', [1 0 0])
errorbar(length(xLocation)+1, jointGrowth, jointConfidence(1,1)-jointGrowth, jointConfidence(1,2)-jointGrowth,'k.');


xticks([xLocation length(xLocation)+1])
xticklabels([metabolitesToPlot {'meanValue'}]);
xtickangle(90)

figure()
hold all

predFlux = allBeta(:,2)/massPerMCell;
confFlux = allConf/massPerMCell;

xLocation = 1:length(predFlux);
bar(xLocation,predFlux)
errorbar(xLocation, predFlux, confFlux(:,1)-predFlux, confFlux(:,2)-predFlux,'k.');
%ylim([-150 50])
xticks(xLocation)
xticklabels(metabolitesToPlot);
xtickangle(90)

%%
labelOut = [metabolitesToPlot {'Growth'}];
fluxOut = [predFlux; jointGrowth];
confOut = [confFlux; jointConfidence];

saveFlux([condition 'SSZ'], labelOut, fluxOut, confOut)
