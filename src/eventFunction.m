function [value,isterminal,direction] = eventFunction(t,y, amounts, haltDir)
    value =(haltDir.*(y-amounts));
    %value(a<0) = 0;
    isterminal = ones(length(amounts),1); %halt for all
    direction = zeros(length(amounts),1);
end
