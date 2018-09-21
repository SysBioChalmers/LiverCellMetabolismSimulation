function [growthRates] = expFit(data, points)

semiLog = data(:,1:2);
semiLog(:,2) = log(semiLog(:,2));

Pfit = findOptima(semiLog, points);

growthRates = cumsum(Pfit);

x = linspace(0, semiLog(end,1));
y = generalizedModel(Pfit, x, points, semiLog(1,:));
        
subplot(2,1,1)    
hold all
title('Growth')
plot(x, exp(y));   
errorbar(data(:,1), data(:,2), data(:,3), 'o')


subplot(2,1,2)
hold all
title('Semi log')
plot(x, y);   
plot(semiLog(:,1),semiLog(:,2),'o')

end


function PfitMin = findOptima(data, points)
n = length(points);
nrOfTests = 15;
ResnormMin = inf;

options = optimset('Display','off');

    for i = 1:nrOfTests
            lb = -0.04 * ones(n+1,1);
            ub =  0.06 * ones(n+1,1);
            P0 = rand(n+1,1) * 0.12 - 0.06;
            model = @(P, x) generalizedModel(P, x, points, data(1,:));
            [Pfit, Resnorm] = lsqcurvefit(model, P0, data(:,1), data(:,2), lb, ub, options);
            if Resnorm<ResnormMin
                ResnormMin = Resnorm;
                PfitMin = Pfit;
            end
    end
end


function model = generalizedModel(slopes, x, breakPoints, startingPoint)
    model = startingPoint(2) + slopes(1) .* x;
    
    for i = 1:length(breakPoints)
        model = model + pf(x-breakPoints(i)).* slopes(i+1);
    end
end


function result = pf(x)
    result = max(x,0);
end

