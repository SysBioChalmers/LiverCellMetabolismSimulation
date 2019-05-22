function composition = getCompositionFromProteomics(file, sequenceProt, sequences)
    proteomicsFile = importdata(file);
    protCount = proteomicsFile.data;
    protCountNames = proteomicsFile.textdata;
    countLookUp = containers.Map(protCountNames, protCount);

    countVector = zeros(size(sequences,1),1);
    hits = 0;
    for i = 1:length(sequenceProt)
        if isKey(countLookUp, sequenceProt{i})
            countVector(i) = countLookUp(sequenceProt{i});
            hits = hits + 1;
        end    
    end
    %hits
    
    countMultiplier = repmat(countVector, 1, size(sequences,2));
    
    composition = countMultiplier .* sequences;
    composition = sum(composition);
    composition = composition/sum(composition);
end

