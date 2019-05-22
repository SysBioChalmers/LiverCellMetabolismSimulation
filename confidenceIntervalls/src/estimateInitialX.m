function x0 = estimateInitialX(xvals, yvals, tvals)
    nrOfTypes = length(unique(xvals(:,2)));
    
    x0 = zeros(nrOfTypes,2);
    
    %estimate biomass by log
    biomassT = xvals(xvals(:,2)==1,1);
    biomassX = yvals(xvals(:,2)==1);

    km = polyfit(biomassT, log(biomassX),1);
    
    x0(1,1) = exp(km(2));
    x0(1,2) = km(1);
    
    %estimate S0 by mean initial concentration
    Sinit = zeros(nrOfTypes-1,1);
    for i = 1:size(Sinit,1)
        match = and(xvals(:,1)==tvals(1), xvals(:,2) == (i+1));
        Sinit(i) = mean(yvals(match));
    end
    
    %estimate f by mean of highest flux
    flux = zeros(nrOfTypes-1,1);
    meanGrowthTerm = x0(1,1)*(exp(x0(1,2)*tvals(end)) - 1)/x0(1,2);
    
    for i = 1:size(flux,1)
        match = xvals(:,2) == (i+1);
        highestT = max(xvals(match,1));
        match = and(xvals(:,1)==highestT, match);
        finalConc = mean(yvals(match));
        flux(i) = (finalConc-Sinit(i))/meanGrowthTerm;
        %finalConc-Sinit(i)
    end    
    
    flux = flux.*(1+0.2.*(rand(length(flux),1)-0.5));
    

    %estimate STD by mean STD from all time points
    std0 = zeros(nrOfTypes,1);
    timePoints = unique(xvals(:,1));
    for i = 1:nrOfTypes
        k = 0;
        for j = 1:length(timePoints)
            match = and(xvals(:,1) == timePoints(j), xvals(:,2) == i);
            if sum(match)>0
                std0(i) = std0(i) + std(yvals(match));
                k = k+1;
            end
        end
        std0(i) = std0(i)/k;
    end
    
    
    x0(2:end,1) = Sinit;
    x0(2:end,2) = flux;
    x0(:,3) = 1.5*std0; %add some margin
end

