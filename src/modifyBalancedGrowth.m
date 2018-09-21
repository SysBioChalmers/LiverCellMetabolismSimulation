function bg = modifyBalancedGrowth(bg, mets)
    gly = findIndex(mets, 'glycine[s]');
    ser = findIndex(mets, 'serine[s]');
    phe = findIndex(mets, 'phenylalanine[s]');
    tyr = findIndex(mets, 'tyrosine[s]');
    ala = findIndex(mets, 'alanine[s]');
    bg(ser) = bg(ser) + bg(gly);
    bg(gly) = 0;
    bg(phe) = bg(phe) + bg(tyr);
    bg(tyr) = 0;
    bg(ala) = 0;
end

