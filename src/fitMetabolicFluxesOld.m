function newFluxValues = fitMetabolicFluxes(model, tspan, fluxMets, fluxValues, halflife, outputMets, init, haltCrit, removalPoints, expData, growthdat)
global proteinPerCell

%Get experimental data
tMax= 100;

options = optimset('Display','iter');
options = optimset(options, 'UseParallel', true);

newFluxValues = fluxValues;

for i = 22
    curData = expData(fluxMets{i});
    tdata = curData(1,:)';
    ydata = curData(2,:)';
    ydata(tdata>tMax) = [];
    tdata(tdata>tMax) = [];
    
    P0 = fluxValues(i,:);

    outputValue = find(ismember(outputMets, fluxMets{i}));
    
    fitFunction = @(P, t) runModel(P, t, fluxValues, i, outputValue, model, tspan, fluxMets, halflife, outputMets, init, haltCrit, removalPoints);
    [Pfit, Resnorm] = lsqcurvefit(fitFunction, P0, tdata, ydata, [], [], options);   
    ydata
    newFluxValues(i,:) = Pfit;    
end


%newFluxValues = lsqnonlin(@(X)minimizationProblem([fluxValues(1:9,:);[X fluxValues(10,2:end)] ;fluxValues(11:end,:)], xdata, ydata, mapping, weight, model, tspan, fluxMets, halflife, outputMets, init, haltCrit, removalPoints),fluxValues(10,1), [], [], options);



'


end


function yResult = runModel(P0, xdata, fluxValues, inputValue, outputValue, model, tspan, fluxMets, halflife, outputMets, init, haltCrit, removalPoints)
    verbose = false;
    fluxValues(inputValue,:) = P0;
    [t, y, yconc] = fullSimulation(model, tspan, fluxMets, fluxValues, halflife, outputMets, init, haltCrit, removalPoints, verbose);
    yResult = interp1(t, yconc(:,outputValue), xdata, 'linear');
end
