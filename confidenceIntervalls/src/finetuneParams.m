function PfitMin = finetuneParams(data, timePoints, mu, init, volumePoints, prior, tolerance)
    options = optimset('Display','off');

    delta = abs(prior * tolerance);
    delta(end) = []; %remove padding 0
    truncPrior = prior;
    truncPrior(end) = [];
    
    lb = -delta;
    ub =  delta;
        
    startGuess = zeros(1,length(delta));
    
    
    
    model = @(P, t) simulateWithBreakPoints(t, timePoints, mu, [truncPrior+P 0], init(1), init(2), volumePoints);
    [Pfit, Resnorm] = lsqcurvefit(model, startGuess, data(:,1), data(:,2), lb, ub, options);
    
    PfitMin = prior + [Pfit 0];
    
end


