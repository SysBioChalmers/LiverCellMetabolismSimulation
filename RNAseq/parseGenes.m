function genes = parseGenes(dataString)
    genes = split(dataString, 'or');
    for i = 1:length(genes)
        genes{i} = cleanString(genes{i});
    end
end

function str = cleanString(str)
    str = strrep(str,'(','');
    str = strrep(str,')','');
    str = strrep(str,' ','');
end