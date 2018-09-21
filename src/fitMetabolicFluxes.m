function newFluxValues = fitMetabolicFluxes(model, tspan, fluxMets, fluxValues, halflife, outputMets, init, haltCrit, volumePoints, expData)
[t, y, yconc, breakPoints, mu] = fullSimulation(model, tspan, fluxMets, fluxValues, halflife, outputMets, init, haltCrit, volumePoints, false);

%Get experimental data
tMax= 140;

newFluxValues = fluxValues;
Xinit = y(1,1);

for i = 1:length(fluxMets)
    fluxMets{i}
    if isKey(expData,fluxMets{i})
        curData = expData(fluxMets{i});
        tdata = curData(1,:)';
        ydata = curData(2,:)';
        ydata(tdata>tMax) = [];
        tdata(tdata>tMax) = [];

        P0 = fluxValues(i,1:(end-1));
        if length(ydata)>(length(P0)-1)
            Pfit = findOptima(P0, [tdata ydata], breakPoints, mu, Xinit, volumePoints);
            newFluxValues(i,1:(end-1)) = Pfit;
        end
    end
end


for i = 1:length(fluxMets)
    fprintf('%s', fluxMets{i})
    for j = 1:size(newFluxValues,2)
        fprintf('\t%f', newFluxValues(i,j));
    end
    fprintf('\n');
end
end

function Pfit = findOptima(P0, data, breakPoints, mu, Xinit, volumePoints)
    lb = zeros(1,length(P0));
    ub = zeros(1,length(P0));
    
    %Can be consumed or produced:
    deltas = diff(data(:,2));
    normDeltas = abs(deltas)/max(abs(deltas));
    deltas(normDeltas>0.01);
    if sign(min(deltas)) == -1
        lb = [];
    end
    if sign(max(deltas)) == 1
        ub = [];
    end

    options = optimset('Display','off');
    model = @(P, t) simulateWithBreakPoints(t, breakPoints, mu, P, Xinit, data(1,2), volumePoints);
    [Pfit, Resnorm] = lsqcurvefit(model, P0, data(:,1), data(:,2), lb, ub, options);

%     hold all
%     plot(data(:,1), data(:,2), 'o');
%     tpoints = linspace(0,180);
%     ydat = simulateWithBreakPoints(tpoints, breakPoints, mu, Pfit, Xinit, data(1,2), volumePoints);
%     plot(tpoints, ydat)
    
end

function yResult = simulateWithBreakPoints(dataT, breakPoints, mu, P, Xinit, Sinit, volumePoints)
    initialX = Xinit;
    initialS = Sinit;
    t = [];
    y = [];
    lastPoint = 0;
    
    for i = 1:length(breakPoints)
        if i == length(breakPoints)
            flux = 0;
        else
            flux = P(i);
        end
        tspan = [lastPoint breakPoints(i)];
        [curt, curConc, ycells] = getMetaboliteCurves(tspan, mu(i), flux, initialX, initialS, volumePoints);
        initialX = ycells(end);
        initialS = curConc(end);        
        t = [t;curt];
        y = [y; curConc];
        lastPoint = breakPoints(i);
    end

    
    %remove duplicate points
    [t, indx] = unique(t); 
    y = y(indx,:);
    yResult = interp1(t, y, dataT, 'linear');
end


function [tVals, yconc, ycells] = getMetaboliteCurves(tspan, mu, flux, initX, initS, volumePoints)
    tVals = [];
    yconc = [];
    ycells = [];
    
    earlierPoints = find(volumePoints(:,1)<=tspan(1));
    earlierPoints(end) = []; %keep the latested concentration
    volumePoints(earlierPoints,:) = [];
    laterPoints = find(volumePoints(:,1)>=tspan(2));
    volumePoints(laterPoints,:) = [];
    
    tpoints = volumePoints(:,1);
    tpoints(1) = tspan(1); %<-start at the T span point
    tpoints(end+1) = tspan(2);

    volumepoints = volumePoints(:,2);
    curInit = [initX, initS];
    
    for i = 1:(length(tpoints)-1)
        curInit(2) = curInit(2) * volumepoints(i);
        if tpoints(i+1)>tspan(2)
            curTspan = [tpoints(i) tspan(2)];
        else
            curTspan = [tpoints(i) tpoints(i+1)];
        end
        
        parameterizedModel = @(t,y) odeModel(t,y,[mu; flux]);
        
        [curT,curY] = ode23s(parameterizedModel, curTspan, curInit);

        %store metabolites as concentration between widhdrawals
        curInit = curY(end,:);
        curInit(2) = curInit(2)/volumepoints(i);
        
        %If halt, report reason
        %StoreResult
        tVals = [tVals; curT];
        yconc = [yconc; curY(:,2)/volumepoints(i)];
        ycells = [ycells; curY(:,1)];
    end
end


function [dy] = odeModel(t,y, rates)
    dy = rates * y(1); %Growth
end

