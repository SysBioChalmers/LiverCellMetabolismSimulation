function fluxes = metTune(data, timePoints, growthRates, initalX, volumePoints, prior, tolerance)
    timePoints = [timePoints'; volumePoints(end,1)];
    init = [initalX; data(1,2)];
    
    fluxes = finetuneParams(data, timePoints, growthRates, init, volumePoints, prior, tolerance);
    representativeFlux = max(abs(fluxes));
    
    if representativeFlux == 0
        precision = 0;
    else
        precision = 2-floor(log10(representativeFlux));
    end
    
    fluxes = round(fluxes, precision);
    
    t = linspace(0,timePoints(end))';
    yold = simulateWithBreakPoints(t, timePoints, growthRates, prior, init(1), init(2), volumePoints);    
    y = simulateWithBreakPoints(t, timePoints, growthRates, fluxes, init(1), init(2), volumePoints);

    hold all
    plot(t, yold, '--');
    plot(t, y);   
    errorbar(data(:,1), data(:,2), data(:,3), 'o')
    ylim([0 inf])
end