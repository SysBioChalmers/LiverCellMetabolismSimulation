clc
aaData = importdata('hepg2-GLN-GLC-AA.txt');
time = aaData.data(1,:);
expData = containers.Map;

for i = 2:length(aaData.data)
    currentMet = aaData.rowheaders{i};
    expData(currentMet) = [time; aaData.data(i,:)];
end

%AA std
aaSTD = importdata('hepg2-GLN-GLC-AA-STD.txt');
for i = 2:length(aaData.data)
    currentMet = aaSTD.rowheaders{i};
    expData(currentMet) = [expData(currentMet); aaSTD.data(i,:)];
end

dict = expData.keys;
for i = 1:length(dict)
   data = expData(dict{i});
   for j = 1:length(data)
        fprintf('%s\t%f\t%f\t%f\n', dict{i}, data(1,j), data(2,j), data(3,j))
   end
end