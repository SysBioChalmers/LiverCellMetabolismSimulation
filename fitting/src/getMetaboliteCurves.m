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
    
    [tVals indx] = unique(tVals);
    ycells = ycells(indx);
    yconc = yconc(indx);    
end


function [dy] = odeModel(t,y, rates)
    dy = rates * y(1); %Growth
end