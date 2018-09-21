function [expData, removalPoints, growth] = loadExpdata(folder, organism, condition)
keySet =   {'Sulfate', 'Riboflavin', 'Folate', 'Pantothenate', 'Niacinamide', 'Choline', 'Phosphate', 'Albumin', 'AGPC', 'Urea', 'Amonia', 'Lactate', 'Glucose', 'Pyruvate', 'Pyroglutamate', 'Alanine', 'Arginine', 'Asparagine', 'Aspartic acid', 'Cystine', 'Glutamic acid', 'Glutamine', 'Glycine', 'Histidine', 'Iso-leucine', 'Leucine', 'Lysine', 'Methionine', 'Ornithine', 'Phenylalanine', 'Serine', 'Taurine', 'Threonine', 'Tryptophan', 'Tyrosine', 'Valine'};
valueSet = {'sulfate[s]', 'riboflavin[s]', 'folate[s]', 'pantothenate[s]', 'nicotinamide[s]', 'choline[s]', 'Pi[s]', 'albumin[s]', 'sn-glycerol-3-PC[s]', 'urea[s]', 'NH3[s]', 'L-lactate[s]', 'glucose[s]', 'pyruvate[s]', '5-oxoproline[s]', 'alanine[s]', 'arginine[s]', 'asparagine[s]', 'aspartate[s]', 'cystine[s]', 'glutamate[s]', 'glutamine[s]', 'glycine[s]', 'histidine[s]', 'isoleucine[s]', 'leucine[s]', 'lysine[s]', 'methionine[s]', 'ornithine[s]', 'phenylalanine[s]', 'serine[s]', 'taurine[s]', 'threonine[s]', 'tryptophan[s]', 'tyrosine[s]', 'valine[s]'};
metTranslation = containers.Map(keySet,valueSet);

pathing = [folder '/' organism '-' condition '-'];
expData = containers.Map;

%Met data

metData = importdata([pathing 'mets.txt']);

for i = 2:length(metData.textdata)
    currentMet = metData.textdata{i,1};  
    currentMet = metTranslation(currentMet);
    currentData = metData.data(i-1,:)';
    if not(expData.isKey(currentMet))
        expData(currentMet) = currentData;
    else
        expData(currentMet) = [expData(currentMet) currentData];
    end
end

addStdValues = true;
if addStdValues
    stdFile = importdata([folder '/standardDeviations.txt']);
    stdMets = stdFile.textdata;
    stdData = stdFile.data;
    for i = 1:length(stdMets)
        curData = expData(stdMets{i});
        curData(3,:) = stdData(i);
        expData(stdMets{i}) = curData;
    end

%Growth data
growthData = importdata([folder '/growthData.txt']);
growthLabels = growthData.textdata(:,1);
growthLabels(1) = []; %Remove header
growthRows = ismember(growthLabels, [organism '-' condition '-']);
growth = growthData.data(growthRows,:);


removalPoints = importdata([folder '/' organism '-' 'volume.txt']); %sample removal points
end

