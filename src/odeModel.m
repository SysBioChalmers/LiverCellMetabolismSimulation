function [dy] = odeModel(t,y, rates, halflife)
    dy = rates * y(1) +  halflife * y; %Growth and chemestry
    
end
