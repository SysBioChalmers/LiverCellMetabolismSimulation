addpath('../../src')
[num,txt,raw] = xlsread('AA_Hepg2_Huh7.xlsx');
filename = 'serialized.txt';

raw(1:2,:) = []; %remove NaN row

hepg2data = false;

%Remove Trailing data
raw(25:end,:) = [];

%Remove NaN cols
removeCols = [];
for i = 2:size(raw,2)
   if isnan(raw{1,i})
       removeCols = [removeCols i];
   end
end
raw(:,removeCols) = [];

if hepg2data
    %Remove Huh7 data
    raw(:,17:end) = [];
else
    %Remove Hepg2 data
    raw(:,2:16) = [];    
end
    
%Merge ID
for i = 2:size(raw,2)
    raw{1,i} = sprintf('%s-%s', num2str(raw{1,i}), num2str(raw{2,i}));
end
raw(2,:) = [];

for i = 1:size(raw,1)
    if strcmp(raw{i,1}, 'iso-leucine')
        raw{i,1} = 'Iso-leucine';
    elseif strcmp(raw{i,1}, 'leucine')
        raw{i,1} = 'Leucine';
    end
end

rawData = raw(4:end,2:end);
rawMets = raw(4:end,1);
rawID = raw(1,2:end);
rawCondition = raw(2,2:end);
rawTime = raw(3,2:end);


fileID = fopen(filename,'w');

%AA time condition sampleId value
fprintf(fileID,'met\ttime\tcondition\tID\tvalue\n');

%Print all mets
for i = 1:length(rawMets)
    curMet = rawMets{i};
    for j = 1:size(rawData,2)
        curTime = num2str(rawTime{j});
        curCondition = num2str(rawCondition{j});
        curID = rawID{j};
        curVal = num2str(rawData{i,j});
        if strcmp(curVal, 'nd')
            curVal = 'NaN';
        end
        %AA time condition sampleId value
        fprintf(fileID,'%s\t%s\t%s\t%s\t%s\n', curMet, curTime, curCondition, curID, curVal);
    end
end

fclose(fileID);