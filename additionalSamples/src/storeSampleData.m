function storeSampleData(raw, samplekey, filename)
    mapRaw = IO(samplekey);
    sampleMap.cols = mapRaw(1,:);
    sampleMap.data = mapRaw(2:end,:);
    sampleMap.ID = cell2nummat(sampleMap.data(:,ismember(sampleMap.cols,'id')));
    rawEntities = raw(:,1);
    rawData = raw(:,2:end);
    rawID = cell2mat(rawData(ismember(rawEntities,'id'),:));
    rawData(ismember(rawEntities,'id'),:) = [];
    rawMets = rawEntities(2:end);
    
    fileID = fopen(filename,'w');
    fprintf(fileID,'%s', makeHeaderEntry(sampleMap));
    
    for i = 1:length(rawMets)
        curMet = rawMets{i};
        for j = 1:size(rawData,2)
            curID = rawID(j);
            curVal = rawData{i,j};
            curStr = makeStringEntry(sampleMap, curMet, curID, curVal);
            fprintf(fileID,'\n%s', curStr);
        end
    end
    
    fclose(fileID);
end

function curStr = makeHeaderEntry(sampleMap)
    %AA sampleId time condition value
    curStr = sprintf('met');
    for i = 1:length(sampleMap.cols)
        curStr = sprintf('%s\t%s', curStr, sampleMap.cols{i});
    end
    curStr = sprintf('%s\tvalue', curStr);
end

function curStr = makeStringEntry(sampleMap, curMet, curID, curVal)
    %AA sampleId time condition value
    curStr = sprintf('%s', curMet);
    for i = 1:length(sampleMap.cols)
        matchingId = (curID == sampleMap.ID);
        curStr = sprintf('%s\t%s', curStr, sampleMap.data{matchingId,i});
    end
    curStr = sprintf('%s\t%s', curStr, curVal);
end