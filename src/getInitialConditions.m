function result = getInitialConditions(data, growthdat, timePoint)
global massPerCell
    growthdat = makeMeanValues(growthdat);
    
    result = zeros(length(data.keys)+1,1);
	result(1) = interp1(growthdat(:,1),growthdat(:,2),timePoint) * massPerCell;

    allKeys = data.keys;
    for i = 2:length(result)
       values = data(allKeys{i-1});
       values = makeMeanValues(values')';
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

function output = makeMeanValues(input)
    allTime = input(:,1);
    allTimepoints = sort(unique(allTime));
    output = zeros(length(allTimepoints),size(input,2));
    for i = 1:length(allTimepoints)
        output(i,1) = allTimepoints(i);
        filter = allTime == allTimepoints(i);
        for j = 2:size(input,2)
            output(i,j) = nanmean(input(filter,j));
        end
    end
end