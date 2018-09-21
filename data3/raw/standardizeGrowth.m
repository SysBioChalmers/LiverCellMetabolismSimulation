clc
growth = importdata('growthData.txt');
time = growth.data(1:2:end,1);
expData = containers.Map;

cellCount = zeros(size(growth.data,1)/2, 4);

for i = 1:size(cellCount,1)
    for j = 2:5
        data = growth.data((i*2-1):(i*2),j);
        cellCount(i,j-1) = mean(data); 
        stdData(i,j-1) = std(data);
    end
end

labelOfCondition = growth.textdata(2:5);

for j = 1:length(labelOfCondition)
    for i = 1:size(cellCount,1)
        fprintf('%s\t%f\t%f\t%f\n', labelOfCondition{j}, time(i), cellCount(i,j), stdData(i,j));
    end
end