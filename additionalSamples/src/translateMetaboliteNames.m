function raw = translateMetaboliteNames(raw, fileName)
    %translate metabolites
    aminoAcids = IO(fileName);
    for i = 1:length(aminoAcids)
        curRow = ismember(raw(:,1), aminoAcids{i,1});
        raw{curRow,1} = aminoAcids{i,2};
    end
end

