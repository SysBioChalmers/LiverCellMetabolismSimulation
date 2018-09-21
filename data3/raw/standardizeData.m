

clc
aaData = importdata('glnAA.txt');
time = aaData.data(:,1);
data = aaData.data;
measurements = length(time)/2;

keySet  = {'Glucose', 'Pyruvate', 'CYS', 'GLU', 'GLN', 'GLY', 'HIS', 'ISO', 'LEU', 'LYS', 'MET', 'PHE', 'SER', 'taurine[s]', 'THR', 'TRP', 'TYR', 'VAL'};
valueSet =  {'Glucose', 'Pyruvate', 'Cystine', 'Glutamic acid', 'Glutamine', 'Glycine', 'Histidine', 'Iso-leucine', 'Leucine', 'Lysine', 'Methionine', 'Phenylalanine', 'Serine', 'Taurine', 'Threonine', 'Tryptophan', 'Tyrosine', 'Valine'};
metTranslation = containers.Map(keySet,valueSet);


														


for i = 2:length(aaData.data)
    for j = 1:measurements
        metName = metTranslation(aaData.textdata{i});
        fprintf('%s\t%f\t%f\t%f\n', metName, time(j), data(j,i), data(measurements+j,i))
    end
end

aaData = importdata('glnPyr.txt');
time = aaData.data(:,1);
expData = containers.Map;
data = aaData.data;
measurements = length(time);

for i = 2:length(aaData.textdata)
    for j = 1:measurements
        metName = metTranslation(aaData.textdata{i});
        fprintf('%s\t%f\t%f\t0\n', metName, time(j), data(j,i))
    end
end