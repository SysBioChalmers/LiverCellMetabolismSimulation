clc
close all
addpath('src')
massPerCell = 426.8 * 10^-12; %g per cell
massPerMCell = massPerCell*10^6;

color = [93 155 211
         215 86 40
         238 178 32]/256;

%leaving out 'Asparagine' 'Butyrate', 'Taurine', 'Citrulline' 
metabolitesToPlot = {'Glucose', 'Lactate', 'Pyruvate', 'Alanine' 'Arginine' 'Aspartic acid' 'Cystine' 'Glutamic acid' 'Glutamine' 'Glycine' 'Histidine' 'Iso-leucine' 'Leucine' 'Lysine' 'Methionine' 'Ornithine' 'Phenylalanine' 'Proline' 'Serine' 'Threonine' 'Tryptophan' 'Tyrosine' 'Valine'};

volData = importdata('volume.txt');
volData(:,1) = volData(:,1)-volData(1,1); %move time to 0


condition = '0';
[dataX, dataY, metlabels] = makeDataStructureNew(condition, metabolitesToPlot);
tvals = unique(dataX(:,1))';

[dataYScaled, scalingfactor] = scaleData(dataX, dataY);

x0 = estimateInitialX(dataX, dataY, tvals);
x0 = estimateInitialX(dataX, dataYScaled, tvals);

modelfun = @(x,y) fitFunction(x, y, tvals, volData(:,2));
%[x,R,J,CovB,MSE,EMI] = nlinfit(dataX, dataY, modelfun, x0, 'ErrorModel', 'combined');
[x,R,J,CovB,MSE,EMI] = nlinfit(dataX, dataYScaled, modelfun, x0);


conf = nlparci(x,R,'covar',CovB);
conf = conf.*[[scalingfactor scalingfactor];[scalingfactor scalingfactor]];

figure()
hold all
predGrowth = x(1,2);
confGrowth = conf(size(x,1)+1,:);
confGrowth = conf(size(x,1)+1,:)/scalingfactor(1);

bar(1,predGrowth)
errorbar(1, predGrowth, confGrowth(1)-predGrowth, confGrowth(2)-predGrowth,'k.');
ylim([0 0.05])

figure()
hold all

predFlux =  x(2:end,2)/massPerMCell;
predFlux = scalingfactor(2:end) .* x(2:end,2)/massPerMCell;

confFlux = conf((size(x,1)+2):end,:)/massPerMCell;
xLocation = 1:length(predFlux);
bar(xLocation,predFlux)
errorbar(xLocation, predFlux, confFlux(:,1)-predFlux, confFlux(:,2)-predFlux,'k.');
ylim([-150 50])
xticks(xLocation)
xticklabels(metlabels(2:end));
xtickangle(90)

figure()
trange = linspace(tvals(1),tvals(end));


for i = 1:length(xLocation)
    subplot(4,6,i)
    hold all
    filter = dataX(:,2)==i;

    xrange = [trange' ones(length(trange), 1)*i];

    [ypred,delta] = nlpredci(modelfun, xrange, x, R, 'Covar', CovB, 'MSE', MSE, 'ErrorModelInfo',EMI, 'SimOpt', 'on', 'Alpha', 0.1);
    lower = ypred - delta;
    upper = ypred + delta;
    ul = [lower upper];
    ul = scalingfactor(i)*[lower upper];
    ypred = scalingfactor(i)*ypred;
    
    fillArea(trange, ul, color(1,:))
    scatter(dataX(filter,1), dataY(filter), 'MarkerFaceColor', color(2,:), 'MarkerEdgeColor', color(2,:),'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 1)
    plot(trange, ypred, 'k-', 'linewidth', 2)

    title(metlabels{i})
    maxval = ylim;
    ylim([0 maxval(2)])    
end


