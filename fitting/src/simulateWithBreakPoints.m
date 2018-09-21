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