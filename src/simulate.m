function [Tout, Yout, Yconc, affectedHaltCrit] = simulate(tspan, init, rates, halflife, haltCrit, volumePoints)
    %Identify which medium removal points to include
    earlierPoints = find(volumePoints(:,1)<=tspan(1));
    earlierPoints(end) = []; %keep the latested concentration
    volumePoints(earlierPoints,:) = [];
    laterPoints = find(volumePoints(:,1)>tspan(2));
    volumePoints(laterPoints,:) = [];
    
    if volumePoints(end,1)<tspan(2)
        volumePoints(end+1,1) = tspan(2);
        volumePoints(end,2) = volumePoints(end-1,2);
    end
       
    tpoints = volumePoints(:,1);
    tpoints(1) = tspan(1); %<-start at the T span point
    volumepoints = volumePoints(:,2);
    
    %Initialize variables
    Tout = [];
    Yout = [];
    Yconc = [];
    affectedHaltCrit = zeros(length(haltCrit),1);
    curInit = init*volumepoints(1); %<-from concentration to  
    curInit(1) = init(1); %Ignore biomass
    
  
    %Set up halt points
    
    
    haltDir = sign(init-haltCrit*volumepoints(1));
     
    for i = 1:(length(tpoints)-1)
        volumeAdjustedHalt = haltCrit*volumepoints(i);
        parametrizedEvent = @(t,y) eventFunction(t,y,volumeAdjustedHalt, haltDir);
        options = odeset('Events',parametrizedEvent);        
        
        curTspan = [tpoints(i) tpoints(i+1)];
        
        parameterizedModel = @(t,y) odeModel(t,y,rates, halflife);
        [curT,curY,ye,ie] = ode23s(parameterizedModel, curTspan, curInit, options);
                
        %If halt, report reason
        %StoreResult
        Tout = [Tout; curT];
        Yout = [Yout; curY];
        Yconc = [Yconc; curY/volumepoints(i)];

        if not(isempty(ie))
           ieT = ie';
           affectedHaltCrit = abs(ieT-volumeAdjustedHalt)< 10^-6;
           break
        else
            %Remove for sampling
            curInit = curY(end,:)';
            deltaVolume = volumepoints(i) - volumepoints(i+1);
            curConcentration = curInit/volumepoints(i);
            curInit = curInit - deltaVolume * curConcentration;
            curInit(1) = curY(end,1)'; %exclude biomass            
        end
        
    end
    
    Yconc(:,1) = 0; %concentration does not apply to biomass
end

