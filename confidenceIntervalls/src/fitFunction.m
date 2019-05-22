function y = fitFunction(p,x,t,volume)
    tNr = x(:,1);
    valueNr = x(:,2);
    y = zeros(size(x,1),1);
    volume = [volume(1); volume];
    
    lookuptable = zeros(size(p,1), length(t));
    lookuptable(1,:) = p(1,1) * exp(p(1,2).*t);
    
    %iterative fit curve between each sample removal points
    lookuptable(2:end,1) = p(2:end,1);
    for i = 2:length(t)
        dt = t(i)-t(i-1);
        growthPart = lookuptable(1,i-1)./p(1,2) * (exp(p(1,2).*dt)-1);
        lookuptable(2:end,i) = lookuptable(2:end,i-1) + p(2:end,2)/volume(i) * growthPart;        
    end
    
    %lookuptable(lookuptable<0) = 0;
    
    
    for i = 1:size(x,1)
        if ismember(tNr(i), t)
            %Use lookup table for speed in the fitting
            y(i) = lookuptable(valueNr(i), tNr(i)==t);
        else
            %Calculate continues values for plotting
            y(i) = intermediatePoints(valueNr(i), tNr(i), lookuptable, p, t, volume);
        end
    end
end

function y = intermediatePoints(valueNr, tNr, lookuptable, p, t, volume)
    tpoints = 1:length(t);
    closestLookup = max(tpoints(tNr>t)); %from below
    dt = tNr-t(closestLookup);
    staringValue = lookuptable(valueNr, closestLookup);
    if valueNr == 1
        y = staringValue * exp(p(1,2).*dt);
    else
        growthPart = lookuptable(1,closestLookup)./p(1,2) * (exp(p(1,2).*dt)-1);
        y = lookuptable(valueNr,closestLookup) + p(valueNr,2)/volume(closestLookup+1) * growthPart;
    end
    %y(y<0) = 0;
end


    