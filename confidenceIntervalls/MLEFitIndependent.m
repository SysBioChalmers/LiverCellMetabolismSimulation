clc
close all
addpath('src')
condition = '22';
%confidencebound = 0.341 * 2; %1 std
confidencebound = 0.95; %2 std

massPerCell = 426.8 * 10^-12; %g per cell
massPerMCell = massPerCell*10^6;

color = [93 155 211
         215 86 40
         238 178 32]/256;

%leaving out  'Asparagine' 'Butyrate', 'Taurine', 'Citrulline' 
metabolitesToPlot = {'Glucose' 'Lactate' 'Pyruvate' 'Alanine' 'Arginine' 'Aspartic acid' 'Cystine' 'Glutamic acid' 'Glutamine' 'Glycine' 'Histidine' 'Iso-leucine' 'Leucine' 'Lysine' 'Methionine' 'Ornithine' 'Phenylalanine' 'Proline' 'Serine' 'Threonine' 'Tryptophan' 'Tyrosine' 'Valine'};
%metabolitesToPlot = {'Glucose' 'Lactate' 'Pyruvate'};

%Load data
volData = importdata('volume.txt');
startTime = volData(1,1);
volData(:,1) = volData(:,1) - startTime; %move time to 0
addConstraints = NaN(length(metabolitesToPlot)+1,3);
[dataX, dataY, metOrder] = makeDataStructureNew(condition, metabolitesToPlot);

%Estimating and constraining STD of cells using only growth data
disp('Estimating STD of growth...')
[dataX, dataY] = makeDataStructureNew(condition, {});
tvals = unique(dataX(:,1))';
x0 = estimateInitialX(dataX, dataY, tvals);
modelfun = @(x,y) fitFunction(x, y, tvals, volData(:,2));
[beta, fval] = mLE(dataX, dataY, modelfun, x0, NaN(1,3));
addConstraints(1,3) = beta(1,3);

%Estimating and constraining STD for each metabolite in turn using the data
disp('Estimating STD of each metabolite...')
for i = 2:length(metOrder)
    [dataX, dataY] = makeDataStructureNew(condition, metOrder(i));
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
    [dataX, dataY, mets] = makeDataStructureNew(condition, {curMet});
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
    subplot(6,5,i)
    hold all
    curBeta = [allGrowthBeta(i,:);allBeta(i,:)];
    plotTrajectory(trange, modelfun, condition, startTime, curBeta, metabolitesToPlot(i), color)
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
jointConfidence = [min(confFlux(:,1)) max(confFlux(:,2))];

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

saveFlux(condition, labelOut, fluxOut, confOut)
