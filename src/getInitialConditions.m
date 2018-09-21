function result = getInitialConditions(data, growthdat, timePoint)
global massPerCell
    result = zeros(length(data.keys)+1,1);
	result(1) = interp1(growthdat(:,1),growthdat(:,2),timePoint) * massPerCell;

    allKeys = data.keys;
    for i = 2:length(result)
       values = data(allKeys{i-1});
       if size(values,2) > 1
           result(i) = interp1(values(1,:),values(2,:),timePoint);
           if isnan(result(i))
               result(i) = 0;
           end
       else
           result(i) = values(2);
       end
    end
end

