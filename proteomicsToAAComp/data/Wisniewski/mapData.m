reference = importdata('../AAComp.txt');
refId = reference.textdata;

fileID = fopen('consensus.txt', 'w');
fid = fopen('allKeys.txt', 'r');
tline = fgetl(fid);
while ischar(tline)
    ids = strsplit(tline,';');
    for i = 1:length(ids)
        tmp = strsplit(ids{i}, '-');
        ids{i} =  tmp{1};
    end
        
    matches = find(ismember(ids,refId));
    if not(isempty(matches))
        fprintf(fileID,'%s\n', ids{matches(end)});
    else
        fprintf(fileID,'%s\n', ids{1});
        disp(ids{1})
    end
    tline = fgetl(fid);
end
fclose(fid);
fclose(fileID);