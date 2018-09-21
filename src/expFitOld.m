function expFit(data)
breakPoints = [0 1 2];

x = linspace(0, data(end,1));
y = zeros(length(x),length(breakPoints));


    for i = 1:length(breakPoints)
        Pfit = findOptima(data, breakPoints(i));
        y(:,i) = generalizedModel(Pfit, x, breakPoints(i), data(1,:));
    end
subplot(1,2,1)
    hold all
    for i = 1:length(breakPoints)
        plot(x, y(:,i));   
    end
    plot(data(:,1),data(:,2),'o')
subplot(1,2,2)    
    hold all
    for i = 1:length(breakPoints)
        plot(x, exp(y(:,i)));   
    end
  
    plot(data(:,1),exp(data(:,2)),'o')
end


function PfitMin = findOptima(data, n)
nrOfTests = 100;
ResnormMin = inf;

options = optimset('Display','off');

    for i = 1:nrOfTests
            lb = [ones(n,1); -0.04 * ones(n+1,1)];
            ub = [data(end,1)*ones(n,1); 0.06 * ones(n+1,1)];
            P0 = [data(end,1)*rand(n,1)*0.5; rand(n+1,1) * 0.12 - 0.06];
            model = @(P, x) generalizedModel(P, x, n, data(1,:));
            [Pfit, Resnorm] = lsqcurvefit(model, P0, data(:,1), data(:,2), lb, ub, options);
            if Resnorm<ResnormMin
                PfitMin = Pfit;
            end
    end
end


function model = generalizedModel(P, x, nrOfPoints, startingPoint)
    breakPoints = cumsum(P(1:nrOfPoints));
    slopes =  P((nrOfPoints+1):end);
    model = startingPoint(2) + slopes(1) .* x;
    
    for i = 1:nrOfPoints
        model = model + pf(x-breakPoints(i)).* slopes(i+1);
    end
end


function result = pf(x)
    result = max(x,0);
end

