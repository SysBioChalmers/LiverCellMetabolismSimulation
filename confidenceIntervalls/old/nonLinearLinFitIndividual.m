clc
close all
addpath('src')
massPerCell = 426.8 * 10^-12; %g per cell
massPerMCell = massPerCell*10^6;

color = [93 155 211
         215 86 40
         238 178 32]/256;

%leaving out 'Butyrate', 'Taurine', 'Citrulline' 
%metabolitesToPlot = {'Glucose', 'Lactate', 'Pyruvate', 'Asparagine', 'Butyrate', 'Taurine', 'Citrulline', 'Alanine' 'Arginine' 'Aspartic acid' 'Cystine' 'Glutamic acid' 'Glutamine' 'Glycine' 'Histidine' 'Iso-leucine' 'Leucine' 'Lysine' 'Methionine' 'Ornithine' 'Phenylalanine' 'Proline' 'Serine' 'Threonine' 'Tryptophan' 'Tyrosine' 'Valine'};
metabolitesToPlot = {'Glucose', 'Lactate', 'Pyruvate', 'Asparagine', 'Alanine' 'Arginine' 'Aspartic acid' 'Cystine' 'Glutamic acid' 'Glutamine' 'Glycine' 'Histidine' 'Iso-leucine' 'Leucine' 'Lysine' 'Methionine' 'Ornithine' 'Phenylalanine' 'Proline' 'Serine' 'Threonine' 'Tryptophan' 'Tyrosine' 'Valine'};

volData = importdata('volume.txt');
volData(:,1) = volData(:,1)-volData(1,1); %move time to 0


%condition = '22mM';
condition = '0';
%leaving out 'Asparagine', 'Butyrate', 'Taurine', 'Citrulline'

trange = linspace(0, 25);

metIndx = 2; %1 = growth 2 = metabolite


predFlux = zeros(length(metabolitesToPlot),1);
confFlux = zeros(length(metabolitesToPlot),2);
    
for i = 1:length(metabolitesToPlot)
    subplot(5,6,i)

    [dataX, dataY, metlabels] = makeDataStructureNew(condition, metabolitesToPlot{i});
    [dataYScaled, scalingfactor] = scaleData(dataX, dataY);
    tvals = unique(dataX(:,1))';
    
    x0 = estimateInitialX(dataX, dataYScaled, tvals);
    modelfun = @(x,y) fitFunction(x, y, tvals, volData(:,2));
    %[x,R,J,CovB,MSE,ErrorModelInfo] = nlinfit(dataX, dataY, modelfun, x0);
    [x,R,J,CovB,MSE,ErrorModelInfo] = nlinfit(dataX, dataYScaled, modelfun, x0);

    conf = nlparci(x,R,'jacobian',J);
    conf = conf.*[[scalingfactor scalingfactor];[scalingfactor scalingfactor]];
    scalingfactor = scalingfactor(metIndx);
    
    predFlux(i) = scalingfactor.* x(metIndx,2)/massPerMCell;
    confFlux(i,:) = conf(2 + metIndx,:)/massPerMCell;
    
    
    hold all
    filter = dataX(:,2)==metIndx;

    xrange = [trange' ones(length(trange), 1)*metIndx];

    [ypred,delta] = nlpredci(modelfun, xrange, x, R, 'Covar', CovB, 'MSE', MSE, 'SimOpt', 'on', 'Alpha', 0.1);

    lower = ypred - delta;
    upper = ypred + delta;
    ul = [lower upper];
    ul = scalingfactor*[lower upper];
    ypred = scalingfactor*ypred;
    
    tvals = unique(dataX(:,1))';
    fillArea(trange, ul, color(1,:))
    scatter(dataX(filter,1), dataY(filter), 'MarkerFaceColor', color(2,:), 'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 1)
    plot(trange, ypred, 'k-', 'linewidth', 2)

    title(metabolitesToPlot{i})
    maxval = ylim;
    ylim([0 maxval(2)])    
end

figure()
hold all
xLocation = 1:length(predFlux);
bar(xLocation,predFlux)
errorbar(xLocation, predFlux, confFlux(:,1)-predFlux, confFlux(:,2)-predFlux,'k.');
ylim([-150 50])
xticks(xLocation)
xticklabels(metabolitesToPlot);
xtickangle(90)


saveFlux(condition, metabolitesToPlot, predFlux, confFlux)
