function plotTrajectory(trange, modelfun, condition, startTime, beta, metName, color)
    hold all

    [dataX, dataY, metlabels] = makeDataStructureNew(condition, metName);
    curMet = findIndex(metlabels, metName);
    filter = dataX(:,2)==curMet;
    xrange = [trange' curMet * ones(length(trange),1)];
    
    ypred = modelfun(beta, xrange);
    delta = beta(curMet,3) * 1.96;
    
    lower = ypred - delta;
    upper = ypred + delta;
    ul = [lower upper];

    fillArea(startTime + trange, ul, color(1,:))
    scatter(startTime + dataX(filter,1), dataY(filter), 'MarkerFaceColor', color(2,:), 'MarkerEdgeColor', color(2,:),'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 1)
    plot(startTime+ trange, ypred, 'k-', 'linewidth', 2)

    title(metName)
    maxval = ylim;
    %ylim([0 maxval(2)])    
end

