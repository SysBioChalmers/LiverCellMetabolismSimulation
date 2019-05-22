networkBg = imread('network.png');
imshow(networkBg);

[groupNames, reactionGroups, rxnStochiometry, coordinates] = importReactionGroups('metaboliteMap.txt');
coordinates = zeros(length(groupNames),2);

for i = 1:length(groupNames)
    disp(groupNames{i})
    coordinates(i,:) = ginput(1);
    text(coordinates(i,1), coordinates(i,2), groupNames{i});
end

coordinates = round(coordinates);

%Store Coordinates
fileID = fopen('metaboliteMap.txt','w');
fprintf(fileID,'GroupName\tRXNS\tStochiometry\tCoordinates');

for i = 1:length(groupNames)
    curGroup = groupNames{i};
    curRxns = strjoin(reactionGroups{i}, ';');
    curStochiometry = strjoin(strsplit(num2str(rxnStochiometry{i})), ';');
    
    curCoordinates = strjoin(strsplit(num2str(coordinates(i,:))), ';');
    fprintf(fileID,'\n%s\t%s\t%s\t%s', curGroup, curRxns, curStochiometry, curCoordinates);
end

fclose(fileID);