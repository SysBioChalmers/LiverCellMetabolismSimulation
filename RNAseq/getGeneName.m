function geneName = getGeneName(geneConvertionList, key)
    indx = ismember(geneConvertionList(:,1), key);
    geneName = unique(geneConvertionList(indx,2));
    geneName = geneName{1};
end

