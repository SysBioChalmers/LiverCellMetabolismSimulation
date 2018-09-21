clc
aaData = importdata('huh7-22mm-AA.txt');
time = aaData.data(1,:);
expData = containers.Map;

for i = 2:length(aaData.data)
    currentMet = aaData.rowheaders{i};
    currentData = [aaData.data(i,:); zeros(1, length(time))];
    expData(currentMet) = [time; currentData];
end

%AA std
aaGlu = importdata('huh7-22mm-GluLact.txt');
time = aaGlu.data(1,:);
for i = 2:size(aaGlu.data,1)
    currentMet = aaGlu.rowheaders{i};
    currentData = [aaGlu.data(i,:); zeros(1, length(time))];
    expData(currentMet) = [time; currentData];
end

fprintf('Name\tTime\tConcentration\tSTD\n')
dict = expData.keys;
for i = 1:length(dict)
   data = expData(dict{i});
   for j = 1:length(data)
        fprintf('%s\t%f\t%f\t%f\n', dict{i}, data(1,j), data(2,j), data(3,j))
   end
end
fprintf('Tryptophan\t%f\t%f\t%f\n', 0, 70, 0)
