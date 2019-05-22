function PfitMin = findOptima(data, timePoints, mu, init, volumePoints)
n = length(mu);
nrOfTests = 5;
ResnormMin = inf;

options = optimset('Display','off');

    for i = 1:nrOfTests
            lb = -1*ones(n,1);
            ub =  ones(n,1);
            P0 = rand(n,1) * 2 - 1;
            model = @(P, t) simulateWithBreakPoints(t, timePoints, mu, P, init(1), init(2), volumePoints);
            [Pfit, Resnorm] = lsqcurvefit(model, P0, data(:,1), data(:,2), lb, ub, options);
            if Resnorm<ResnormMin
                ResnormMin = Resnorm;
                PfitMin = Pfit;
            end
    end
end



