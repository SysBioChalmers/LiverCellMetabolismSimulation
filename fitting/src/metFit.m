function fluxes = metFit(data, timePoints, growthRates, initalX, volumePoints)
    timePoints = [timePoints'; volumePoints(end,1)];
    init = [initalX; data(1,2)];

    fluxes = findOptima(data, timePoints, growthRates, init, volumePoints);

    t = linspace(0,timePoints(end))';
    y = simulateWithBreakPoints(t, timePoints, growthRates, fluxes, init(1), init(2), volumePoints);

    hold all
    plot(t, y);   
    errorbar(data(:,1), data(:,2), data(:,3), 'o')
    ylim([0 inf])
end