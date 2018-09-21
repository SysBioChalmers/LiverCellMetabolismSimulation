function [t, y, yconc, breakPoints, mu] = fullSimulation(model, tspan, fluxMets, fluxValues, halflife, outputMets, init, haltCrit, volumePoints, verbose)
    t = [];
    y = [];
    yconc = [];
    breakPoints =[];
    mu = [];

    for i = 1:size(fluxValues, 2)
        model = bindFBA(model, fluxMets, fluxValues(:,i)/1000);
        
        rates = runFBA(model, outputMets, verbose);
        
        
        rates(2:end) = rates(2:end) * 1000;
        
        [curt, curY, curYconc, affectedHaltCrit] = simulate(tspan, init, rates, halflife, haltCrit, volumePoints);
        
                
        %Move halt crit to 0 if above. Only halt once for crit.
        if haltCrit(affectedHaltCrit==1) == 0
            haltCrit(affectedHaltCrit==1) = -1000;
        else
            haltCrit(affectedHaltCrit==1) = 0;
        end
        
        if verbose
            disp(outputMets(affectedHaltCrit == 1))
        end
        
        breakPoints = [breakPoints; curt(end)];
        mu = [mu;rates(1)];
        
        if length(curYconc)>1
            t = [t;curt];
            y = [y; curY];
            yconc = [yconc; curYconc];

            init = curYconc(end,:)';
            init(1) = curY(end,1);
            tspan(1) = curt(end);
            
            if tspan(1)>=tspan(2)
                break
            end
        else
            break 
        end
    end

%remove duplicate points
[t, indx] = unique(t); 
y = y(indx,:);
yconc = yconc(indx,:);    
end

