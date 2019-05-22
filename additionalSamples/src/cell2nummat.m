function nummat = cell2nummat(cellarray)
    nummat = cellfun(@str2num, cellarray);
end

