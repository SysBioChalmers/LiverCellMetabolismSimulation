close all
cellLine = {
    'hepg2'
    'hepg2'
    'hepg2' 
    'huh7'
    'huh7'    
    'huh7'
    };

condition = {
    '0'
    '6'
    '22' 
    '0'
    '6'
    '22' 
    };

metabolites = {
    'Alanine'      
    'Arginine'     
    'Aspartic acid'
    'Cystine'      
    'Glutamic acid'
    'Glutamine'    
    'Glycine'      
    'Histidine'    
    'Iso-leucine'  
    'Leucine'      
    'Lysine'       
    'Methionine'   
    'Ornithine'    
    'Phenylalanine' 
    'Serine'       
    'Threonine'    
    'Tyrosine'     
    'Valine'};

metConcentration = zeros(length(metabolites),length(condition));


metaboliteMap = IO('../data/metaboliteMap.tsv');
metTranslation = containers.Map(metaboliteMap(:,1),metaboliteMap(:,2));






for i = 1:length(condition)
    [expData] = loadExpdata('../data', cellLine{i}, condition{i});
    
    for j = 1:length(metabolites)
        curMet = metTranslation(metabolites{j});
        if isKey(expData, curMet)
            curData = expData(curMet);
            expDataFirst = curData(2,1);
            metConcentration(j,i) = expDataFirst;
        end
    end
end
standardDev = nanstd(metConcentration')'/1000; %unit mM


%Cells
cellLine = {
    'hepg2'
    'hepg2'
    'hepg2'
    };

condition = {
    '0'
    '6'
    '22'};

timePoints = [0 54.5 71.5];
allSTD = zeros(length(condition), length(timePoints));

for i = 1:length(condition)
    [expData, removalPoints, growthData] = loadExpdata('../data', cellLine{i}, condition{i});
    for j = 1:length(timePoints)
        filter = growthData(:,1) == timePoints(j);
        curData = growthData(filter,2);
        allSTD(i,j) = nanstd(curData);
    end
end

metabolites = [metabolites; {'growth'}];
standardDev = [standardDev; mean(allSTD(:))/10^6]; %unit per million

%Glucose
timePoints = [54.25 71.5];

allSTD = zeros(length(condition), length(timePoints));
for i = 1:length(condition)
    expData = loadExpdata('../data', cellLine{i}, condition{i});
    glucoseData = expData('glucose[s]');
    for j = 1:length(timePoints)
        filter = glucoseData(1,:) == timePoints(j);
        curData = glucoseData(2,filter);
        allSTD(i,j) = nanstd(curData);
    end
end
allSTD = allSTD(:);
allSTD(allSTD==0) = [];

metabolites = [metabolites; {'Glucose'}];
standardDev = [standardDev; mean(allSTD)/1000]; %unit mM

%Lactate
allSTD = zeros(length(condition), length(timePoints));
for i = 1:length(condition)
    expData = loadExpdata('../data', cellLine{i}, condition{i});
    lactateData = expData('L-lactate[s]');
    for j = 1:length(timePoints)
        filter = lactateData(1,:) == timePoints(j);
        curData = lactateData(2,filter);
        allSTD(i,j) = nanstd(curData);
    end
end
allSTD = allSTD(:);
allSTD(allSTD==0) = [];

metabolites = [metabolites; {'Lactate'}];
standardDev = [standardDev; mean(allSTD)/1000]; %unit mM

%Pyruvate
allSTD = zeros(length(condition), length(timePoints));
for i = 1:length(condition)
    expData = loadExpdata('../data', cellLine{i}, condition{i});
    pyruvateData = expData('pyruvate[s]');
    for j = 1:length(timePoints)
        filter = pyruvateData(1,:) == timePoints(j);
        curData = pyruvateData(2,filter);
        allSTD(i,j) = nanstd(curData);
    end
end
allSTD = allSTD(:);
allSTD(allSTD==0) = [];

metabolites = [metabolites; {'Pyruvate'}];
standardDev = [standardDev; mean(allSTD)/1000]; %unit mM


filename = 'StandardDeviations.tsv';

fileID = fopen(filename,'w');
fprintf(fileID,'Met\tSTD');

for i = 1:length(metabolites)
    fprintf(fileID,'\n%s\t%2.5f\t%2.5f', metabolites{i}, standardDev(i));
end

fclose(fileID);
