clc
close all
addpath('src')
massPerCell = 426.8 * 10^-12; %g per cell
massPerMCell = massPerCell*10^6;

color = [93 155 211
         215 86 40
         238 178 32]/256;

%leaving out 'Asparagine', 'Butyrate', 'Taurine', 'Citrulline' 
metabolitesToPlot = {'Glucose', 'Lactate', 'Pyruvate', 'Alanine' 'Arginine' 'Aspartic acid' 'Cystine' 'Glutamic acid' 'Glutamine' 'Glycine' 'Histidine' 'Iso-leucine' 'Leucine' 'Lysine' 'Methionine' 'Ornithine' 'Phenylalanine' 'Proline' 'Serine' 'Threonine' 'Tryptophan' 'Tyrosine' 'Valine'};
%metabolitesToPlot = {'Glutamic acid'};

volData = importdata('volume.txt');
volData(:,1) = volData(:,1)-volData(1,1); %move time to 0

condition = '0';
[dataX, dataY, metlabels] = makeDataStructureNew(condition, metabolitesToPlot);
tvals = unique(dataX(:,1))';

priorSTD = loadPriorSTD('StandardDeviations.tsv', metlabels);

x0 = estimateInitialX(dataX, dataY, tvals);

std0 = zeros(length(metlabels),1);
for i = 1:length(metlabels)
    if priorSTD(i) == 0
        filter = and(dataX(:,1) == 0, dataX(:,2) == i);
        std0(i) = std(dataY(filter));
    else
        std0(i) = priorSTD(i);
    end
end

priorSTD(1) = 0;
priorSTD(2) = 0;

modelfun = @(x,y) fitFunction(x, y, tvals, volData(:,2));
[beta, fval] = mLE(dataX, dataY, modelfun, x0, std0);

conf = mLEConfidence(dataX, dataY, modelfun, beta);

%%
figure()
trange = linspace(tvals(1),tvals(end));
for i = 1:length(metlabels)
    subplot(5,6,i)
    hold all
    filter = dataX(:,2)==i;

    xrange = [trange' ones(length(trange), 1)*i];

    %[ypred,delta] = nlpredci(modelfun, xrange, beta, R, 'Covar', CovB, 'MSE', MSE, 'SimOpt', 'on', 'Alpha', 0.1);
    ypred = modelfun(beta, xrange);
    delta = beta(i,3) * 1.96;
    
    lower = ypred - delta;
    upper = ypred + delta;
    ul = [lower upper];

    fillArea(trange, ul, color(1,:))
    scatter(dataX(filter,1), dataY(filter), 'MarkerFaceColor', color(2,:), 'MarkerEdgeColor', color(2,:),'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 1)
    plot(trange, ypred, 'k-', 'linewidth', 2)

    title(metlabels{i})
    maxval = ylim;
    ylim([0 maxval(2)])    
end

%%

%Find the std of that fits 95% confidence interval
xrange = [trange' ones(length(trange), 1)*1];
ypred = modelfun(beta, xrange);
