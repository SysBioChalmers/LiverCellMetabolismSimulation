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

% function yResult = getMetaboliteCurves(t, timePoints, mu, flux, curInit)
% tVals = [];
% yconc = [];
%     for i = 1:length(mu)
%         curTspan = [timePoints(i) timePoints(i+1)];
%         currentTpoints = t;
%         currentTpoints(currentTpoints<curTspan(1)) = [];
%         currentTpoints(currentTpoints>curTspan(2)) = [];
%         curTspan = [curTspan(1);currentTpoints;curTspan(2)];
%         curTspan =unique(curTspan);
%         parameterizedModel = @(t,y) odeModel(t,y,[mu(i); flux(i)]);
%         [curT,curY] = ode23s(parameterizedModel, curTspan', curInit);
%         curInit = curY(end,:);
%         tVals = [tVals;curT];
%         yconc = [yconc;curY(:,2)];
%     end
% 
% yResult = zeros(length(t),1);
%     for i = 1:length(t)
%         match = find(tVals == t(i));
%         yResult(i) = yconc(match(1));
%     end
%        
% end

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

