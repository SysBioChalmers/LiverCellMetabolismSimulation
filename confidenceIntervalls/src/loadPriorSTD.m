function priorSTD = loadPriorSTD(fileName, metlabels)
    priorSTD = zeros(length(metlabels),1);
    stdData = IO(fileName);
    mets = stdData(2:end,1);
    vals = cell2nummat(stdData(2:end,2));
    
    for i = 1:length(metlabels)
       curMet = findIndex(mets, metlabels{i}); 
       if not(isempty(curMet))
            priorSTD(i) = vals(curMet);
       end
    end
end

